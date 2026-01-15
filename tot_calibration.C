#include <map>
#include <vector>
#include <iosfwd>
#include <istream>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

#include <TFile.h>
#include <TTree.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TF1.h>
#include <TLatex.h>
#include <TProfile.h>
#include <TGraph.h>
#include <TSpline.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TParameter.h>
#include <TError.h>

#include "common.C"

float energy_fraction_cut = 0.3;
long n_events = 1000000;
// n_events = 100;


float tot_min_cut = 0000;

std::vector<int> run_numbers = {330};
std::vector<float> energies  = {3.8};


void tot_calibration() {
    gErrorIgnoreLevel = kWarning;

    run_numbers.clear();
    energies.clear();
    for (int r = 319; r <= 338; r++) {
        run_numbers.push_back(r);
        energies.push_back(1.6 + 0.2 * (r - 319));
    }

    gStyle->SetOptStat(0);
    auto mapping = read_mapping("eeemcal_desy_dec2025_mapping.csv");

    float adc_calib = 0;
    TFile *adc_calib_file = TFile::Open("output/adc_to_gev_calibration.root", "READ");
    if (adc_calib_file && !adc_calib_file->IsZombie()) {
        TParameter<float>* adc_calib_param = (TParameter<float>*)adc_calib_file->Get("mean_adc_to_gev_calibration");
        if (adc_calib_param) {
            adc_calib = adc_calib_param->GetVal();
            std::cout << "Loaded ADC to GeV calibration: " << adc_calib << std::endl;
        }
    }
    if (adc_calib == 0) {
        std::cerr << "Error: ADC to GeV calibration not found!" << std::endl;
        return;
    }
    adc_calib_file->Close();

    TH1* gain_factor = nullptr;
    TH1* crystal_gain_factor;
    TFile* gain_file = TFile::Open("output/gain_factors.root");
    if (gain_file && !gain_file->IsZombie()) {
        gain_factor = (TH1*)gain_file->Get("gain_factors");
        crystal_gain_factor = (TH1*)gain_file->Get("crystal_factor");
        std::cout << "Loaded gain factors from file." << std::endl;
    }
    if (!gain_factor) {
        gain_factor = new TH1F("gain_factors", "Gain Factor", 400, 0, 400);
        for (int i = 1; i <= 400; i++) {
            gain_factor->SetBinContent(i, 1.0);
        }
    }
    if (!crystal_gain_factor) {
        crystal_gain_factor = new TH1F("crystal_factor", "Crystal Gain", 25, 0, 25);
        for (int i = 1; i <= 25; i++) {
            crystal_gain_factor->SetBinContent(i, 1);
        }
    }

    std::vector<TH1*> missing_energies;
    std::vector<TH1*> tot_distributions;
    std::vector<TH1*> tot_energies;
    std::vector<TH2*> cog_distributions;
    std::vector<TH1*> tot_bin_distribution;
    std::vector<TH2*> tot_bin_vs_value;
    std::vector<std::vector<TH1*>> per_sipm_tot;
    std::vector<std::vector<std::vector<TH2*>>> tot_covariance;

    TH2 *total_distribution =new TH2F("run%d_tot_energy", "ToT Energy;Energy (GeV);ToT Signal", 100, 0, 5, 1000, 0, 60000);
    TH2 *total_distribution_invt =new TH2F("run%d_tot_energy_invt", "Energy as a function of ToT;ToT Signal;Energy (GeV);", 1000, 0, 60000, 100, 0, 5);
    


    bool use_widths = false;
    TFile *input_file = TFile::Open("output/tot_widths.root", "READ");
    if (!input_file || input_file->IsZombie()) {
        std::cerr << "Error: Could not open input file 'output/tot_width.root'" << std::endl;
    } else {
        use_widths = true;
    }

    for (int run = 0; run < run_numbers.size(); run++) {
        TH1 *peak_widths = nullptr;
        if (use_widths) {
            peak_widths = (TH1*)input_file->Get(Form("run%d_tot_widths", run_numbers[run]));
            if (!peak_widths) {
                std::cerr << "Warning: Could not find histogram 'run" << run_numbers[run] << "_tot_widths' in input file." << std::endl;
                use_widths = false;
            }
            else {
                std::cout << "Using ToT widths for run " << run_numbers[run] << std::endl;
            }
        }


        int run_number = run_numbers[run];
        float energy = energies[run];
        // Process data file
        TFile* root_file = TFile::Open(Form("/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run%03d.root", run_number));
        TTree* tree = (TTree*)root_file->Get("events");
        uint32_t adc[576][20];
        uint32_t tot[576][20];
        int tot_events = 0;
        tree->SetBranchAddress("adc", &adc);
        tree->SetBranchAddress("tot", &tot);
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("adc", 1);
        tree->SetBranchStatus("tot", 1);

        
        missing_energies.push_back(new TH1F(Form("run%d_missing_energy", run_number), Form("Run %d ADC portion of energy", run_number), 100, 0, 6));
        tot_distributions.push_back(new TH1F(Form("run%d_tot_distribution", run_number), Form("Run %d ToT;ToT;Events", run_number), 1024, 0, 16 * 4096));
        tot_energies.push_back(new TH2F(Form("run%d_tot_energy", run_number), Form("Run %d ToT Energy;Energy (GeV);(ToT Signal)", run_number), 100, 0, 5, 1000, 0, 60000));
        cog_distributions.push_back(new TH2F(Form("run%d_cog_distribution", run_number), Form("Run %d Center of Gravity Distribution;X (# Crystals);Y (# Crystals)", run_number), 100, -0.5, 4.5, 100, -0.5, 4.5));
        tot_bin_distribution.push_back(new TH1F(Form("run%d_tot_bin_distribution", run_number), Form("Run %d ToT Bin Distribution;ToT Bin;Events", run_number), 20, 0, 20));
        tot_bin_vs_value.push_back(new TH2F(Form("run%d_tot_bin_vs_value", run_number), Form("Run %d ToT Bin vs ToT Value;ToT Value;ToT Bin", run_number), 20, 0, 20, 1024, 0, 4096));

        per_sipm_tot.push_back(std::vector<TH1*>());
        tot_covariance.push_back(std::vector<std::vector<TH2*>>());
        for (int i = 0; i < sipms_to_use; i++) {
            per_sipm_tot[run].push_back(new TH1F(Form("run%d_sipm%d_tot_distribution", run_number, i), Form("Run %d SiPM %d ToT;ToT;Events", run_number, i), 1024, 0, 4096));
            tot_covariance[run].push_back(std::vector<TH2*>());
            for (int j = 0; j < sipms_to_use; j++) {
                tot_covariance[run][i].push_back(new TH2F(Form("run%d_sipm%d_sipm%d_tot_covariance", run_number, i, j),
                                                        Form("Run %d SiPM %d vs SiPM %d ToT;SiPM %d ToT;SiPM %d ToT", run_number, i, j, i, j),
                                                        1024, 0, 4096, 1024, 0, 4096));
            }
        }

        Long64_t nentries = tree->GetEntries();
        if (n_events < nentries) {
            nentries = n_events;
        }
        std::cout << "Run " << run_number << std::endl;
        std::cout << "Processing " << nentries << " events" << std::endl;
        int complete = 0;
        for (Long64_t entry = 0; entry < nentries; ++entry) {
            if (entry * 25 / nentries > complete) {
                complete = entry * 25 / nentries;
                print_progress(complete);
            }

            tree->GetEntry(entry);
            float central_signal = 0.0f;
            float central_nine_signal = 0.0f;
            float total_signal = 0.0f;
            bool is_tot_event = false;

            float signals[25];
            for (int crystal = 0; crystal < 25; crystal++) {
                if (crystal == 9) {
                    continue;
                }
                if (crystal == 12) {
                    continue;
                }
                float crystal_signal = 0.0f;
                for (int sipm = 0; sipm < sipms_to_use; sipm++) {
                    int channel = mapping[crystal][sipm];
                    float gain = gain_factor->GetBinContent(crystal * 16 + sipm + 1);
                    float channel_signal = calculate_signal(adc[channel], gain);
                    crystal_signal += channel_signal;
                    if (!is_tot_event && is_tot(tot[channel])) {
                        is_tot_event = true;
                    }
                }
                signals[crystal] = crystal_signal * crystal_gain_factor->GetBinContent(crystal + 1);
            }
            if (is_tot_event) {
                tot_events++;
                continue;
            }

            // Get the ToT
            float center_signal = 0;
            std::vector<float> sipm_signals(16, 0);
            int channels_used = 0;
            bool all_tot = true;
            for (int sipm = 0; sipm < sipms_to_use; sipm++) {
                float signal = 0;
                int channel = mapping[12][sipm];
                float gain = gain_factor->GetBinContent(12 * 16 + sipm + 1);
                all_tot &= calculate_signal(adc[channel], tot[channel], gain, signal);
                if (!all_tot) {
                    break;
                }
                if (use_widths) {
                    float channel_mean = peak_widths->GetBinContent(sipm + 1);
                    float channel_sigma = peak_widths->GetBinError(sipm + 1);
                    // std::cout << "signal: " << signal << ", mean: " << channel_mean << ", sigma: " << channel_sigma << std::endl;
                    if (signal < channel_mean - (2 * channel_sigma) || signal > channel_mean + (2 * channel_sigma)) {
                        all_tot = false;
                        break;
                        continue;
                    }
                }
                sipm_signals[sipm] = signal;
                    // tot_bin_distribution[run]->Fill(timebin);
                    // tot_bin_vs_value[run]->Fill(timebin, signal);
                per_sipm_tot[run][sipm]->Fill(signal);  
                center_signal += signal * gain;
                channels_used++;
            }
            if (!all_tot || channels_used != sipms_to_use) {
                continue;
            }
            // normalize the center signal to 16 channels used
            if (channels_used > 0) {
                center_signal *= ((float)sipms_to_use / channels_used);
            }

            if (center_signal < tot_min_cut) {
                center_signal = 0;
            }
            // Fill the covariance matrices
            for (int i = 0; i < sipms_to_use; i++) {
                for (int j = 0; j < sipms_to_use; j++) {
                    tot_covariance[run][i][j]->Fill(sipm_signals[i], sipm_signals[j]);
                }
            }


            float x_cog, y_cog;
            bool keep = calculate_cog(cog_distributions[run], signals);
            keep &= (!is_tot_event);    // Get rid of events with a tot in a non-central channel
            if (!keep) {
                continue;
            }

            tot_distributions[run]->Fill(center_signal);

            // Sum up the total energy found in the non-central crystals
            float non_central_energy = 0;
            for (int crystal = 0; crystal < 25; crystal++) {
                if (crystal == 12) {continue;}
                non_central_energy += signals[crystal];
            }
            // Calibrate to GeV
            non_central_energy /= adc_calib;
            missing_energies[run]->Fill(non_central_energy);
            tot_energies[run]->Fill(energy - non_central_energy, center_signal);
            total_distribution->Fill(energy - non_central_energy, center_signal);
            total_distribution_invt->Fill(center_signal, energy - non_central_energy);
        }
        // std::cout << std::endl;
        std::cout << "\rTotal ToT events skipped: " << tot_events << std::endl;
    }

    
    if (use_widths) {
        input_file->Close();
    }


    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    text->SetTextFont(42);

    std::vector<std::vector<double>> peak_mean;
    std::vector<std::vector<double>> peak_sigma;

    TCanvas* canvas = new TCanvas("tot_calib", "", 800, 600);
    canvas->SaveAs("output/tot_calib.pdf(");
    for (int run = 0; run < run_numbers.size(); run++) {
        peak_mean.push_back(std::vector<double>());
        peak_sigma.push_back(std::vector<double>());
        canvas->Clear();
        canvas->Divide(3, 2);
        canvas->cd(1);
        missing_energies[run]->Draw();
        text->SetTextAlign(31);
        text->DrawLatex(0.85, 0.85, Form("Mean: %.2f", missing_energies[run]->GetMean()));
        text->DrawLatex(0.85, 0.80, Form("ADC percentage: %.2f", 100 * missing_energies[run]->GetMean() / energies[run]));

        canvas->cd(2);
        tot_distributions[run]->Draw();
        canvas->cd(3);
        tot_energies[run]->Draw("colz");
        canvas->cd(4);
        cog_distributions[run]->Draw("colz");

        auto pad = canvas->cd(5);
        pad->Divide(4, 4, 0, 0);
        for (int sipm = 0; sipm < sipms_to_use; sipm++) {
            pad->cd(sipm + 1);
            // pad->SetMargin(0, 0, 0, 0);
            per_sipm_tot[run][sipm]->SetTitle("");
            float mean = per_sipm_tot[run][sipm]->GetBinCenter(per_sipm_tot[run][sipm]->GetMaximumBin());
            float rms = 500;
            float fit_min = mean - 1.5 * rms;
            float fit_max = mean + 1.5 * rms;
            auto fit = new TF1("tot_fit", "gaus", fit_min, fit_max);
            fit->SetParameter(0, per_sipm_tot[run][sipm]->GetBinContent(per_sipm_tot[run][sipm]->GetMaximumBin()));
            fit->SetParameter(1, mean);
            fit->SetParameter(2, rms);
            auto result = per_sipm_tot[run][sipm]->Fit(fit, "RQS");
            // check if the fit failed
            if (result->Status() != 0 || fit->GetNDF() < 20) {
                std::cerr << "Warning: Fit failed for run " << run_numbers[run] << " SiPM " << sipm << std::endl;
                std::cout << "       Using default fit parameters." << std::endl;
                delete fit;
                fit = new TF1("tot_fit", "gaus", 2000, 3000);
                fit->SetParameter(0, 200);
                fit->SetParameter(1, 2200);
                fit->SetParameter(2, 500);
                per_sipm_tot[run][sipm]->Fit(fit, "RQS");
            }
            per_sipm_tot[run][sipm]->Draw("hist");
            fit->Draw("same");
            peak_mean[run].push_back(fit->GetParameter(1));
            peak_sigma[run].push_back(fit->GetParameter(2));
        }

        canvas->cd(6);
        tot_bin_vs_value[run]->Draw("colz");

        canvas->SaveAs("output/tot_calib.pdf");

        bool do_covariance = false;
        if (do_covariance) {
            canvas->Clear();
            canvas->Divide(16, 16, 0, 0);
            for (int i = 0; i < sipms_to_use; i++) {
                for (int j = 0; j < sipms_to_use; j++) {
                    if (j <= i) {
                        canvas->cd(i * 16 + j + 1);
                        tot_covariance[run][i][j]->Draw("colz");
                    }
                }
            }
            canvas->SaveAs("output/tot_calib.pdf");
        }
    }
    canvas->Clear();
    total_distribution->Draw("colz");
    canvas->SaveAs("output/tot_calib.pdf");

    float scale_factor = 16.0f / sipms_to_use;

    TF1 *range_one_fit = new TF1("range_one", "pol1", 35000 / scale_factor, 50000 / scale_factor);
    total_distribution_invt->Fit(range_one_fit, "RQ");
    range_one_fit->SetLineColor(kRed);

    TF1 *range_two_fit = new TF1("range_two", "pol1", 20000 / scale_factor, 26000 / scale_factor);
    total_distribution_invt->Fit(range_two_fit, "RQ");
    range_two_fit->SetLineColor(kGreen + 2);

    TF1 *total_fit = new TF1("total_fit", "pol2", 20000 / scale_factor, 55000 / scale_factor);
    total_distribution_invt->Fit(total_fit, "RQ");
    total_fit->SetLineColor(kMagenta);


    total_distribution_invt->Draw("colz");
    range_one_fit->Draw("same");
    range_two_fit->Draw("same");
    total_fit->Draw("same");
    canvas->SaveAs("output/tot_calib.pdf)");

    TFile *tot_calib_file = TFile::Open("output/tot_calibration_values.root", "RECREATE");
    TParameter<float>* tot_c0 = new TParameter<float>("tot_c0", total_fit->GetParameter(0));
    TParameter<float>* tot_c1 = new TParameter<float>("tot_c1", total_fit->GetParameter(1));
    TParameter<float>* tot_c2 = new TParameter<float>("tot_c2", total_fit->GetParameter(2));
    TParameter<float>* tot_a0 = new TParameter<float>("tot_a0", range_two_fit->GetParameter(0));
    TParameter<float>* tot_a1 = new TParameter<float>("tot_a1", range_two_fit->GetParameter(1));
    tot_c0->Write();
    tot_c1->Write();
    tot_c2->Write();
    tot_a0->Write();
    tot_a1->Write();
    tot_calib_file->Close();

    if (!use_widths) {
        TFile *tot_widths = TFile::Open("output/tot_widths.root", "RECREATE");
        for (int run = 0; run < run_numbers.size(); run++) {
            TH1* tot_width = new TH1F(Form("run%d_tot_widths", run_numbers[run]), Form("Run %d ToT Widths;SiPM;Width (ADC counts)", run_numbers[run]), 16, 0, 16);
            for (int sipm = 0; sipm < sipms_to_use; sipm++) {
                float mean = peak_mean[run][sipm];
                float sigma = peak_sigma[run][sipm];
                tot_width->SetBinContent(sipm + 1, mean);
                tot_width->SetBinError(sipm + 1, sigma);
            }
            tot_width->Write();
        }
    }
}