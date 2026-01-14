#include <map>
#include <vector>
#include <iostream>
#include <fstream>
#include <istream>
#include <iosfwd>
#include <sstream>

#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2.h>
#include <TH2F.h>
#include <TLatex.h>
#include <TLine.h>
#include <TPad.h>
#include <TStyle.h>
#include <TTree.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TGraphErrors.h>
#include <TParameter.h>

#include "common.C"

float energy_fraction_cut = 0.3;
long n_events = 1000000;
// n_events = 1000;
float tot_min_cut = 0;


void energy_resolution() {

    gStyle->SetOptStat(0);
    
    std::vector<int> run_numbers = {316, 321, 326, 331, 336};
    std::vector<float> energies  = {1.0, 2.0, 3.0, 4.0, 5.0};
    run_numbers.clear();
    energies.clear();
    for (int r = 316; r <= 338; r++) {
        // if (r == 318 || r == 319) {
        //     continue;
        // }
        run_numbers.push_back(r);
        energies.push_back(1. + 0.2 * (r - 316));
        // if (r == 316) {
        //     r += 5;
        // }
    }

    float adc_calib = 0;
    bool values_set = true;
    TFile *adc_calib_file = TFile::Open("output/adc_to_gev_calibration.root", "READ");
    if (adc_calib_file && !adc_calib_file->IsZombie()) {
        TParameter<float>* adc_calib_param = (TParameter<float>*)adc_calib_file->Get("mean_adc_to_gev_calibration");
        if (adc_calib_param) {
            adc_calib = adc_calib_param->GetVal();
            std::cout << "Loaded ADC to GeV calibration: " << adc_calib << std::endl;
        } else {
            values_set = false;
        }
    } else {
        values_set = false;
    }
    if (!values_set) {
        std::cerr << "Error: ADC to GeV calibration not found!" << std::endl;
        return;
    }
    adc_calib_file->Close();

    values_set = true;
    float c0, c1, c2, a0, a1 = 0;
    TFile *tot_calib_file = TFile::Open("output/tot_calibration_values.root", "READ");
    if (tot_calib_file && !tot_calib_file->IsZombie()) {
        TParameter<float>* c0_param = (TParameter<float>*)tot_calib_file->Get("tot_c0");
        TParameter<float>* c1_param = (TParameter<float>*)tot_calib_file->Get("tot_c1");
        TParameter<float>* c2_param = (TParameter<float>*)tot_calib_file->Get("tot_c2");
        TParameter<float>* a0_param = (TParameter<float>*)tot_calib_file->Get("tot_a0");
        TParameter<float>* a1_param = (TParameter<float>*)tot_calib_file->Get("tot_a1");
        if (c0_param && c1_param && c2_param && a0_param && a1_param) {
            c0 = c0_param->GetVal();
            c1 = c1_param->GetVal();
            c2 = c2_param->GetVal();
            a0 = a0_param->GetVal();
            a1 = a1_param->GetVal();
            std::cout << "Loaded ToT calibration: " << c0 << ", " << c1 << ", " << c2 << ", " << a0 << ", " << a1 << std::endl;
        } else {

            values_set = false;
        }
    } else {
        values_set = false;
    }
    if (!values_set) {
        std::cerr << "Error: ToT calibration not found!" << std::endl;
        return;
    }
    tot_calib_file->Close();


    int central_crystal_index = 12;
    int center_nine_indexes[8] = {7, 8, 9, 11, 13, 17, 18, 19};
    int remaining_indexes[16] = {0, 1, 2, 3, 4, 5, 6, 10, 11, 15, 16, 20, 21, 22, 23, 24};

    std::vector<TH1*> central_crystal_energy_vec;
    std::vector<TH1*> central_nine_energy_vec;
    std::vector<TH1*> total_energy_vec;
    std::vector<TH2*> cog_distribution_vec;
    std::vector<TH1*> tot_fraction_vec;

    std::vector<std::vector<TH1*>> crystal_energy;
    std::vector<std::vector<TH1*>> crystal_energy_shares;

    bool use_widths = false;
    TFile *input_file = TFile::Open("output/tot_widths.root", "READ");
    if (!input_file || input_file->IsZombie()) {
        std::cerr << "Error: Could not open input file 'output/tot_width.root'" << std::endl;
    } else {
        use_widths = true;
    }

    // Set up all the histograms we need
    for (int run = 0; run < run_numbers.size(); run++) {
        int run_number = run_numbers[run];
        central_crystal_energy_vec.push_back(new TH1F(Form("run%d_central_crystal_energy", run_number), Form("Run %d Central Crystal Energy;Energy (ADC);Events", run_number), 400, 0, 8));
        central_nine_energy_vec.push_back(new TH1F(Form("run%d_central_nine_energy", run_number), Form("Run %d Central 3x3 Energy;Energy (ADC);Events", run_number), 400, 0, 8));
        total_energy_vec.push_back(new TH1F(Form("run%d_total_energy", run_number), Form("Run %d Total Energy;Energy (ADC);Events", run_number), 400, 0, 8));
        cog_distribution_vec.push_back(new TH2F(Form("run%d_cog_distribution", run_number), Form("Run %d Center of Gravity Distribution;X (# Crystals);Y (# Crystals)", run_number), 100, -0.5, 4.5, 100, -0.5, 4.5));
        tot_fraction_vec.push_back(new TH1F(Form("run%d_tot_fraction", run_number), Form("Run %d ToT Fraction;ToT Fraction;Events", run_number), 17, 0, 17));

        crystal_energy.push_back(std::vector<TH1*>());
        crystal_energy_shares.push_back(std::vector<TH1*>());
        for (int i = 0; i < 25; i++) { 
            crystal_energy[run].push_back(new TH1F(Form("run%d_crystal_%02d_energy", run_number, i), Form("Run %d Crystal %02d Energy;Energy (ADC);Events", run_number, i), 500, 0, 8));
            crystal_energy_shares[run].push_back(new TH1F(Form("run%d_crystal_%02d_energy_share", run_number, i), Form("Run %d Crystal %02d Energy Share;Energy Share;Events", run_number, i), 10000, 0, 1));
        }
    }

    auto mapping = read_mapping("eeemcal_desy_dec2025_mapping.csv");

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

    for (int run = 0; run < run_numbers.size(); run++) {
        int events_displayed = 0;
        int events_used = 0;
        // Process data file
        int run_number = run_numbers[run];
        float energy = energies[run];
        char file_path[256];
        sprintf(file_path, "/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run%03d.root", run_number);
        TFile* root_file = TFile::Open(file_path);
        TTree* tree = (TTree*)root_file->Get("events");
        uint32_t adc[576][20];
        uint32_t tot[576][20];
        int tot_events = 0;
        tree->SetBranchAddress("adc", &adc);
        tree->SetBranchAddress("tot", &tot);
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("adc", 1);
        tree->SetBranchStatus("tot", 1);

        TH1 *peak_widths = nullptr;
        // if (use_widths) {
        //     peak_widths = (TH1*)input_file->Get(Form("run%d_tot_widths", run_numbers[run]));
        //     if (!peak_widths) {
        //         std::cerr << "Warning: Could not find histogram 'run" << run_numbers[run] << "_tot_widths' in input file." << std::endl;
        //     }
        //     else {
        //         std::cout << "Using ToT widths for run " << run_numbers[run] << std::endl;
        //     }
        // }

        Long64_t nentries = tree->GetEntries();
        if (n_events < nentries) {
            nentries = n_events;
        }
        std::cout << "Run " << run_number << " at " << energy << " GeV:" << std::endl;
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
                for (int sipm = 0; sipm < 16; sipm++) {
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

            // Get the center
            float center_signal = 0;
            bool center_is_tot = false;
            bool mixed_event = false;
            int num_tot = 0;
            for (int sipm = 0; sipm < 16; sipm++) {
                float signal = 0;
                int channel = mapping[12][sipm];
                float gain = gain_factor->GetBinContent(12 * 16 + sipm + 1);

                bool tot_sig = calculate_signal(adc[channel], tot[channel], gain, signal);
                if (tot_sig) {
                    num_tot++;
                }
                if (tot_sig && peak_widths) {
                    float channel_mean = peak_widths->GetBinContent(sipm + 1);
                    float channel_sigma = peak_widths->GetBinError(sipm + 1);
                    // std::cout << "signal: " << signal << ", mean: " << channel_mean << ", sigma: " << channel_sigma << std::endl;
                    if (signal < channel_mean - (2 * channel_sigma) || signal > channel_mean + (2 * channel_sigma)) {
                        mixed_event = true;
                        // break;
                        continue;
                    }
                }

                center_is_tot |= tot_sig;
                center_signal += signal;
            }

            tot_fraction_vec[run]->Fill(num_tot);

            if (false && num_tot == 8 && events_displayed < 10) {
                events_displayed++;
                std::vector<TH2*> sipm_waveforms;
                for (int sipm = 0; sipm < 16; sipm++) {
                    sipm_waveforms.push_back(new TH2F(Form("run%d_event%d_sipm%d_waveform", run_number, entry, sipm),
                                                     Form("Run %d Event %d SiPM %d Waveform;Time Bin;ADC Signal", run_number, entry, sipm),
                                                     20, 0, 20, 1024, 0, 1024));
                    int channel = mapping[12][sipm];
                    for (int t = 0; t < 20; t++) {
                        sipm_waveforms[sipm]->Fill(t, adc[channel][t]);
                    }
                    for (int t = 0; t < 10; t++) {
                        if (is_tot(tot[channel])) {
                            sipm_waveforms[sipm]->Fill(7., 200.); 
                        } else {
                            // sipm_waveforms[sipm]->Fill(0., 0.);
                        }
                    }
                    for (int i = 0; i < 20; i++) {
                        sipm_waveforms[sipm]->Fill(0., 0.);
                    }
                }
                TCanvas* c1 = new TCanvas("c1", "c1", 1200, 800);
                c1->Divide(4, 4);
                for (int sipm = 0; sipm < 16; sipm++) {
                    c1->cd(sipm + 1);
                    sipm_waveforms[sipm]->Draw("COLZ");
                }
                c1->SaveAs(Form("output/run%d_event%d_center_tot_waveforms.png", run_number, entry));
                c1->Close();
                for (int sipm = 0; sipm < 16; sipm++) {
                    delete sipm_waveforms[sipm];
                }
            }


            if (num_tot > 0 && num_tot < 16) {
                mixed_event = true;
            }

            if (energy >= 1.9 && num_tot == 0) {
                continue;
            }

            if (mixed_event) {
                continue;
            }

            for (int i = 0; i < 25; i++) {
                if (i == 12) {
                    continue;
                }
                signals[i] /= adc_calib;
            }
            
            if (center_is_tot) {
                // Convert the ToT to energy - these are _incredibly_ rough
                float center_energy = 0;
                float poly_energy = c0 + c1 * center_signal + c2 * center_signal * center_signal;
                float lin_energy = a0 + a1 * center_signal;
                if (lin_energy > poly_energy) {
                    center_energy = lin_energy;
                } else {
                    center_energy = poly_energy;
                }
                signals[12] = center_energy;
            } else {
                // std::cout << "center is not tot!" << std::endl;
                signals[12] = center_signal * crystal_gain_factor->GetBinContent(13);
                signals[12] /= adc_calib;
                // std::cout << "center energy: " << signals[12] << std::endl;
            }

            float x_cog, y_cog;
            bool keep = calculate_cog(cog_distribution_vec[run], signals);
            keep &= (!is_tot_event);
            // if (!keep) {
            //     std::cout << "skipping?" << std::endl;
            //     continue;
            // }


            // Populate energies
            central_signal = signals[central_crystal_index];
            for (int i = 0; i < 8; i++) {
                int idx = center_nine_indexes[i];
                central_nine_signal += signals[idx];
            }
            for (int i = 0; i < 16; i++) {
                int idx = remaining_indexes[i];
                total_signal += signals[idx];
            }
            if (signals[12] / total_signal < energy_fraction_cut) { // Skip events where the majority of the energy share is not in the central crystal;
                continue;
            }

            events_used++;

            central_nine_signal += central_signal;
            total_signal += central_nine_signal;
            central_crystal_energy_vec[run]->Fill(central_signal);
            central_nine_energy_vec[run]->Fill(central_nine_signal);
            total_energy_vec[run]->Fill(total_signal);
            // std::cout << total_signal << std::endl;

            // Energy share fill
            for (int i = 0; i < 25; i++) {
                crystal_energy[run][i]->Fill(signals[i]);
                crystal_energy_shares[run][i]->Fill(signals[i] / total_signal);
            }
            
        }
        std::cout << "\nUsed " << events_used << " events" << std::endl;
        fit_peak(central_crystal_energy_vec[run]);
        fit_peak(central_nine_energy_vec[run]);
        fit_peak(total_energy_vec[run]);
    }


    int crystal_mapping[25] = {
        4, 9, 14, 19, 24,
        3, 8, 13, 18, 23,
        2, 7, 12, 17, 22,
        1, 6, 11, 16, 21,
        0, 5, 10, 15, 20
    };

    TH1 *dummy = new TH1F("energy_resolution_center", "Energy Resolution;Energy (GeV);#sigma(E)/E", 1, 0, 6);
    TGraphErrors *res_center = new TGraphErrors(run_numbers.size()-2);
    TGraphErrors *res_middle = new TGraphErrors(run_numbers.size()-2);
    TGraphErrors *res_all    = new TGraphErrors(run_numbers.size()-2);
    res_center->SetMarkerColor(kRed);
    res_center->SetLineColor(kRed);
    res_center->SetLineStyle(kDashed);
    res_center->SetMarkerStyle(21);
    res_center->SetMarkerSize(2);

    res_middle->SetMarkerColor(kBlue);
    res_middle->SetLineColor(kBlue);
    res_middle->SetLineStyle(kDashed);
    res_middle->SetMarkerStyle(22);
    res_middle->SetMarkerSize(2);

    res_all->SetMarkerColor(kMagenta);
    res_all->SetLineColor(kMagenta);
    res_all->SetLineStyle(kDashed);
    res_all->SetMarkerStyle(23);
    res_all->SetMarkerSize(2);

    int point = 0;
    for (int run = 0; run < run_numbers.size(); run++) {
        if (run_numbers[run] == 318 || run_numbers[run] == 319) {
            continue;
        }
        // std::cout << "a run " << run_numbers[run] << std::endl;
        
        // Center crystal resolution with uncertainty
        float sigma = central_crystal_energy_vec[run]->GetFunction("final_fit")->GetParameter(2);
        float mu = central_crystal_energy_vec[run]->GetFunction("final_fit")->GetParameter(1);
        float sigma_err = central_crystal_energy_vec[run]->GetFunction("final_fit")->GetParError(2);
        float mu_err = central_crystal_energy_vec[run]->GetFunction("final_fit")->GetParError(1);
        float res = sigma / mu;
        float res_err = res * sqrt((sigma_err/sigma)*(sigma_err/sigma) + (mu_err/mu)*(mu_err/mu));
        res_center->SetPoint(point, energies[run], 100.0 * res);
        res_center->SetPointError(point, 0, 100.0 * res_err);
        
        // Central 9 resolution with uncertainty
        sigma = central_nine_energy_vec[run]->GetFunction("final_fit")->GetParameter(2);
        mu = central_nine_energy_vec[run]->GetFunction("final_fit")->GetParameter(1);
        sigma_err = central_nine_energy_vec[run]->GetFunction("final_fit")->GetParError(2);
        mu_err = central_nine_energy_vec[run]->GetFunction("final_fit")->GetParError(1);
        res = sigma / mu;
        res_err = res * sqrt((sigma_err/sigma)*(sigma_err/sigma) + (mu_err/mu)*(mu_err/mu));
        res_middle->SetPoint(point, energies[run], 100.0 * res);
        res_middle->SetPointError(point, 0, 100.0 * res_err);
        
        // Total energy resolution with uncertainty
        sigma = total_energy_vec[run]->GetFunction("final_fit")->GetParameter(2);
        mu = total_energy_vec[run]->GetFunction("final_fit")->GetParameter(1);
        sigma_err = total_energy_vec[run]->GetFunction("final_fit")->GetParError(2);
        mu_err = total_energy_vec[run]->GetFunction("final_fit")->GetParError(1);
        res = sigma / mu;
        res_err = res * sqrt((sigma_err/sigma)*(sigma_err/sigma) + (mu_err/mu)*(mu_err/mu));
        res_all->SetPoint(point, energies[run], 100.0 * res);
        res_all->SetPointError(point, 0, 100.0 * res_err);
        
        point++;
    }

    TCanvas* canvas = new TCanvas("gain_matching", "", 800, 600);
    dummy->Draw();
    dummy->SetMinimum(0);
    dummy->SetMaximum(10);
    res_center->Draw("same lp");
    res_middle->Draw("same lp");
    res_all->Draw("same lp");

    TLegend *l = new TLegend(0.6, 0.6, 0.89, 0.89);
    l->SetLineWidth(0);
    l->AddEntry(res_center, "Center Crystal");
    l->AddEntry(res_middle, "Central 9 Crystals");
    l->AddEntry(res_all, "All Crystals");
    l->Draw();

    TF1 *fit_func = new TF1("energy_res_fit", "sqrt([0]*[0] + [1]*[1]/x + [2]*[2]/(x*x))", 1, 5.4);
    fit_func->SetParameter(0, 2);
    fit_func->SetParameter(1, 2);
    fit_func->SetParameter(2, 4);
    res_all->Fit(fit_func, "R");
    fit_func->Draw("same");

    auto text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.03);
    text->SetTextColor(kBlack);
    text->DrawLatex(0.15, 0.30, "#sigma/E = A #oplus B/#sqrt{E} #oplus C/E");
    text->DrawLatex(0.15, 0.25, Form("A: %.2f#pm%.2f %%", fit_func->GetParameter(0), fit_func->GetParError(0)));
    text->DrawLatex(0.15, 0.20, Form("B: %.2f#pm%.2f %%", fit_func->GetParameter(1), fit_func->GetParError(1)));
    text->DrawLatex(0.15, 0.15, Form("C: %.2f#pm%.2f %%", fit_func->GetParameter(2), fit_func->GetParError(2)));

    
    canvas->SaveAs("output/energy_resolution.pdf(");
    canvas->Clear();

    // res_center->Clear();
    // res_middle->Clear();
    // res_all->Clear();

    // for (int run = 0; run < run_numbers.size(); run++) {
    //     float res = central_crystal_energy_vec[run]->GetFunction("final_fit")->GetParameter(2) / central_crystal_energy_vec[run]->GetFunction("final_fit")->GetParameter(1);
    //     res /= energies[run];
    //     res_center->SetPoint(run, energies[run], 100.0 * res);
    //     res = central_nine_energy_vec[run]->GetFunction("final_fit")->GetParameter(2) / central_nine_energy_vec[run]->GetFunction("final_fit")->GetParameter(1);
    //     res /= energies[run];
    //     res_middle->SetPoint(run, energies[run], 100.0 * res);
    //     res = total_energy_vec[run]->GetFunction("final_fit")->GetParameter(2) / total_energy_vec[run]->GetFunction("final_fit")->GetParameter(1);
    //     res /= energies[run];
    //     res_all->SetPoint(run, energies[run], 100.0 * res);
    // }

    // dummy->GetYaxis()->SetTitle("#sigma(E) / E (%)");

    // dummy->Draw();
    // dummy->SetMinimum(0);
    // dummy->SetMaximum(10);
    // res_center->Draw("same lp");
    // res_middle->Draw("same lp");
    // res_all->Draw("same lp");

    // TLegend *l2 = new TLegend(0.6, 0.6, 0.89, 0.89);
    // l2->SetLineWidth(0);
    // l2->AddEntry(res_center, "Center Crystal");
    // l2->AddEntry(res_middle, "Central 9 Crystals");
    // l2->AddEntry(res_all, "All Crystals");
    // l2->Draw();

    // canvas->SaveAs("output/energy_resolution.pdf(");
    // canvas->Clear();



    TH1 *reconstructed_energy = new TH1F("reconstructed_energy", "Reconstructed Energy;Beam Energy (GeV);Reconstructed Energy (GeV);Events", 1, 0, 6);
    TGraph *reconstructed_energy_graph = new TGraph();
    reconstructed_energy_graph->SetMarkerColor(kBlack);
    reconstructed_energy_graph->SetLineColor(kBlack);
    reconstructed_energy_graph->SetMarkerStyle(21);
    reconstructed_energy_graph->SetMarkerSize(2);
    point = 0;
    for (int i = 0; i < run_numbers.size(); i++) {
        // std::cout << "b run " << run_numbers[i] << std::endl;
        if (run_numbers[i] == 318 || run_numbers[i] == 319) {
            continue;
        }

        float beam_energy = energies[i];
        float peak_energy = total_energy_vec[i]->GetFunction("final_fit")->GetParameter(1);
        reconstructed_energy_graph->SetPoint(point, beam_energy, peak_energy);
        point++;
    }
    canvas->SetGrid();
    reconstructed_energy->SetMinimum(0);
    reconstructed_energy->SetMaximum(6);
    reconstructed_energy->Draw("g");
    reconstructed_energy_graph->Draw("same p");

    canvas->SaveAs("output/energy_resolution.pdf");
    canvas->SetGrid(0, 0);



    for (int run = 0; run < run_numbers.size(); run++) {
        // std::cout << "c run " << run_numbers[run] << std::endl;

        canvas->Clear();
        int run_number = run_numbers[run];
        float energy = energies[run];

        canvas->SetRightMargin(0.05);
        central_crystal_energy_vec[run]->Draw("HIST e");
        auto fit = central_crystal_energy_vec[run]->GetFunction("final_fit");
        fit->Draw("same");
        draw_text(fit, run_number, energy);
        canvas->SaveAs("output/energy_resolution.pdf");

        central_nine_energy_vec[run]->Draw("HIST e");
        fit = central_nine_energy_vec[run]->GetFunction("final_fit");
        fit->Draw("same");
        draw_text(fit, run_number, energy);
        canvas->SaveAs("output/energy_resolution.pdf");

        total_energy_vec[run]->Draw("HIST e");
        fit = total_energy_vec[run]->GetFunction("final_fit");
        fit->Draw("same");
        draw_text(fit, run_number, energy);
        canvas->SaveAs("output/energy_resolution.pdf");

        canvas->Clear();
        tot_fraction_vec[run]->Draw("HIST e");
        gPad->SetLogy();
        canvas->SaveAs("output/energy_resolution.pdf");
        gPad->SetLogy(0);


        canvas->SetRightMargin(0.1);
        cog_distribution_vec[run]->Draw("COLZ");
        gPad->SetLogz();

        std::vector<TLine*> lines;
        for (int i = 1; i < 5; i++) {
            float loc = i - 0.5f;
            TLine* line_x = new TLine(loc, -0.5, loc, 4.5);
            line_x->SetLineStyle(2);
            line_x->Draw();
            TLine* line_y = new TLine(-0.5, loc, 4.5, loc);
            line_y->SetLineStyle(2);
            line_y->Draw();
        }
        TLine* line_x = new TLine(2, -0.5, 2, 4.5);
        line_x->SetLineStyle(2);
        line_x->SetLineColor(kRed);
        line_x->Draw();

        TLine* line_y = new TLine(-0.5, 2, 4.5, 2);
        line_y->SetLineStyle(2);
        line_y->SetLineColor(kRed);
        line_y->Draw();

        TEllipse* circle = new TEllipse(center_x, center_y, sigma_x, sigma_y);
        circle->SetLineColor(kRed);
        circle->SetLineWidth(2);
        circle->SetFillStyle(0);
        circle->Draw();

        canvas->SaveAs("output/energy_resolution.pdf");

        canvas->Clear();
        canvas->Divide(5, 5);
        for (int i = 0; i < 25; i++) {
            canvas->cd(i + 1);
            crystal_energy[run][crystal_mapping[i]]->Draw("HIST e");
            // gPad->SetLogx();
        }

        canvas->SaveAs("output/energy_resolution.pdf");
        

        std::vector<float> gev_calib;
        canvas->Clear();
        canvas->Divide(5, 5);
        float mean_calib = 0;
        for (int i = 0; i < 25; i++) {
            canvas->cd(i + 1);
            crystal_energy_shares[run][crystal_mapping[i]]->Draw("HIST e");
            crystal_energy_shares[run][crystal_mapping[i]]->GetXaxis()->SetRangeUser(0.001, 1);
            gPad->SetLogx();
            float signal_for_1gev = crystal_energy[run][crystal_mapping[i]]->GetMean() / crystal_energy_shares[run][crystal_mapping[i]]->GetMean();
            float x_coord = 0.4;
            if (i == 12) {
                x_coord = 0.15;
            }
            TLatex *text = new TLatex();
            text->SetNDC();
            text->SetTextSize(0.04);
            if (crystal_mapping[i] == 9) {
                continue;
            }
        }
        canvas->SaveAs("output/energy_resolution.pdf");
        mean_calib /= 24;   // Since we are excluding crystal 9

    }
    canvas->Clear();
    canvas->SaveAs("output/energy_resolution.pdf)");

    // float mean_x = cog_distribution->GetMean(1);
    // float mean_y = cog_distribution->GetMean(2);
    // float sigma_x = cog_distribution->GetStdDev(1);
    // float sigma_y = cog_distribution->GetStdDev(2);

    // std::cout << "Mean X: " << mean_x << " +/- " << sigma_x << std::endl;
    // std::cout << "Mean Y: " << mean_y << " +/- " << sigma_y << std::endl;

    // std::cout << "Total TOT events: " << tot_events << " out of " << nentries << std::endl;

}

