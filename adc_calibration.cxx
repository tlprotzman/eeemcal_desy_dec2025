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
#include <TProfile.h>
#include <TParameter.h>
#include <TError.h>

#include "common.C"

float energy_fraction_cut = 0.0;
long n_events = 1000000;
// n_events = 1000;

int run_number = 411;
float beam_energy = 1;


void adc_calibration() {
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    float energy = 1;

    int central_crystal_index = 12;
    int center_nine_indexes[8] = {7, 8, 9, 11, 13, 17, 18, 19};
    int remaining_indexes[16] = {0, 1, 2, 3, 4, 5, 6, 10, 11, 15, 16, 20, 21, 22, 23, 24};
    int common_mode_channel[25] = {224, 296, 188, 386, 152,
                                   512, 548, 242, 134, 404,
                                   206, 368, 332, 458,   8,
                                   440, 116, 350,  62,  98,
                                   422, 170, 260, 314, 278};
    /*
    a - 8
    b - 26
    c - 62
    d - 44
    */


    TH1* central_crystal_energy = new TH1F("central_crystal_energy", "Central Crystal Energy;Energy (`);Events", 500, 0, 75000);
    TH1* central_nine_energy    = new TH1F("central_nine_energy", "Central 3x3 Energy;Energy (ADC);Events", 500, 0, 75000);
    TH1* total_energy           = new TH1F("total_energy", "Total Energy;Energy (ADC);Events", 500, 0, 75000);
    TH2* cog_distribution       = new TH2F("cog_distribution", "Center of Gravity Distribution;X (# Crystals);Y (# Crystals)", 100, -0.5, 4.5, 100, -0.5, 4.5);
    TH1 *toa_distribution        = new TH1F("toa_distribution", "ToA Distribution;ToA;Events", 1024, 0, 1024);
    TH1 *toa_sample = new TH1F("toa_sample", "ToA Sample Distribution;Sample;Events", 20, 0, 20);
    
    TH1 *max_sample_index = new TH1F("max_sample_index", "Max Sample Index Distribution;Sample Index;Events", 20, 0, 20);
    TH1 *second_max_sample_index = new TH1F("second_max_sample_index", "Second Max Sample Index Distribution;Sample Index;Events", 20, 0, 20);
    TH1 *third_max_sample_index = new TH1F("third_max_sample_index", "Third Max Sample Index Distribution;Sample Index;Events", 20, 0, 20);
    
    std::vector<TH2*> E_vs_toa;
    for (int sipm = 0; sipm < sipms_to_use; sipm++) {
        E_vs_toa.push_back(new TH2F(Form("E_vs_toa_sipm_%d", sipm), Form("SiPM %d Energy vs ToA;ToA;Energy (ADC)", sipm), 1024, 0, 1024, 500, 0, 3000));
    }   
    
    std::vector<std::vector<TH2*>> toa_correlations;
    for (int i = 0; i < sipms_to_use; i++) {
        toa_correlations.push_back(std::vector<TH2*>());
        for (int j = 0; j < sipms_to_use; j++) {
            toa_correlations[i].push_back(new TH2F(Form("toa_correlation_sipm_%02d_sipm_%02d", i, j),
                                                  Form("SiPM %02d vs SiPM %02d ToA;SiPM %02d ToA;SiPM %02d ToA", i, j, i, j),
                                                  1024, 0, 1024, 1024, 0, 1024));
        }
    }

    std::vector<std::vector<TH1*>> sipm_energy;
    std::vector<std::vector<TH2*>> sipm_waveform;
    std::vector<TH1*> crystal_energy;
    std::vector<TH1*> crystal_energy_shares;
    std::vector<TH2*> crystal_common_mode;
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < sipms_to_use; j++) {
            sipm_energy.push_back(std::vector<TH1*>());
            sipm_energy[i].push_back(new TH1F(Form("crystal_%02d_sipm_%02d_energy", i, j), Form("Crystal %02d SiPM %02d Energy;Energy (ADC);Events", i, j), 500, 0, 4000));
            sipm_waveform.push_back(std::vector<TH2*>());
            sipm_waveform[i].push_back(new TH2F(Form("crystal_%02d_sipm_%02d_waveform", i, j), Form("Crystal %02d SiPM %02d Waveform;Sample;ADC", i, j), 20, 0, 20, 1024, 0, 1024));
        }   
        crystal_energy.push_back(new TH1F(Form("crystal_%02d_energy", i), Form("Crystal %02d Energy;Energy (ADC);Events", i), 500, 0, 40000));
        crystal_energy_shares.push_back(new TH1F(Form("crystal_%02d_energy_share", i), Form("Crystal %02d Energy Share;Energy Share;Events", i), 10000, 0, 1));
        crystal_common_mode.push_back(new TH2F(Form("crystal_%02d_common_mode", i), Form("Crystal %02d Common Mode;Sample;ADC", i), 20, 0, 20, 1024, 0, 1024));
    }

    TH2* pedestals = new TH2F("pedestals", "Pedestals;Channel;Pedestal (ADC)", 72*8, 0, 72*8, 1024, 0, 1024);

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


    // Process data file
    TFile* root_file = TFile::Open(Form("/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run%03d.root", run_number));
    TTree* tree = (TTree*)root_file->Get("events");
    uint32_t adc[576][20];
    uint32_t tot[576][20];
    uint32_t toa[576][20];
    int tot_events = 0;
    tree->SetBranchAddress("adc", &adc);
    tree->SetBranchAddress("tot", &tot);
    tree->SetBranchAddress("toa", &toa);
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("adc", 1);
    tree->SetBranchStatus("tot", 1);
    tree->SetBranchStatus("toa", 1);

    Long64_t nentries = tree->GetEntries();
    if (n_events < nentries) {
        nentries = n_events;
    }
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
        uint32_t central_toa = 0;
        int index = mapping[12][0]; // Central crystal, first SiPM for event ToA selection
        bool outside_toa_range = false;
        for (int i = 0; i < 20; i++) {
            if (toa[index][i]) {
                central_toa = toa[index][i];
                if (central_toa < 200 || central_toa > 800) {
                    outside_toa_range = true;
                    break;
                }
            }
        }
        if (outside_toa_range) {
            // continue;
        }

        float signals[25];
        for (int crystal = 0; crystal < 25; crystal++) {
            float mean_toa = 0;
            uint32_t toa_used = 0;
            if (crystal == 9) {
                continue;
            }
            float crystal_signal = 0.0f;
            for (int i = 0; i < 20; i++) {
                int channel = common_mode_channel[crystal];
                crystal_common_mode[crystal]->Fill(i, adc[channel][i]);
            }

            for (int sipm = 0; sipm < sipms_to_use; sipm++) {
                int channel = mapping[crystal][sipm];
                float gain = gain_factor->GetBinContent(crystal * 16 + sipm + 1);
                float channel_signal = calculate_signal(adc[channel], gain);
                // float channel_signal = calculate_signal_v6(adc[channel], adc[common_mode_channel[crystal]], gain);
                sipm_energy[crystal][sipm]->Fill(channel_signal);
                int i1, i2, i3 = -1;
                float pedestal = (adc[channel][0] + adc[channel][1] + adc[channel][2]) / 3.0f;
                float max1 = 0.0f, max2 = 0.0f, max3 = 0.0f;
                
                for (int i = 0; i < 20; ++i) {
                    sipm_waveform[crystal][sipm]->Fill(i, adc[channel][i]);
                    // if (i < 5) continue; // Skip pedestal samples
                    // if (i >= 10) continue;
                    float sample = adc[channel][i] - pedestal;
                    if (sample > max1) {
                        i3 = i2;
                        i2 = i1;
                        i1 = i;
                        max3 = max2;
                        max2 = max1;
                        max1 = sample;
                    } else if (sample > max2) {
                        i3 = i2;
                        i2 = i;
                        max3 = max2;
                        max2 = sample;
                    } else if (sample > max3) {
                        i3 = i;
                        max3 = sample;
                    }
                }
                max_sample_index->Fill(i1);
                second_max_sample_index->Fill(i2);
                third_max_sample_index->Fill(i3);




                // pedestals->Fill(channel, (adc[channel][0] + adc[channel][1] + adc[channel][2]) / 3.0f);
                pedestals->Fill(channel,(adc[channel][0]));
                pedestals->Fill(channel,(adc[channel][1]));
                pedestals->Fill(channel,(adc[channel][2]));
                crystal_signal += channel_signal;
                uint32_t this_toa = get_toa(toa[channel]);
                if (this_toa >= 0) {
                    if (crystal == 12) {
                        toa_distribution->Fill(this_toa);
                        int timebin = -1;
                        for (int sample = 0; sample < 20; sample++) {
                            if (toa[channel][sample]) {
                                timebin = sample;
                                break;
                            }
                        }
                        toa_sample->Fill(timebin);
                        E_vs_toa[sipm]->Fill(this_toa, channel_signal);
                        for (int other_sipm = 0; other_sipm < sipms_to_use; other_sipm++) {
                            if (other_sipm == sipm) {
                                continue;
                            }
                            int other_channel = mapping[crystal][other_sipm];
                            uint32_t other_toa = get_toa(toa[other_channel]);
                            if (other_toa < 0) {
                                continue;
                            }
                            toa_correlations[sipm][other_sipm]->Fill(this_toa, other_toa);
                        }
                    }
                    mean_toa += this_toa;
                    toa_used++;
                }
                if (!is_tot_event && is_tot(tot[channel])) {
                    is_tot_event = true;
                }
            }
            signals[crystal] = crystal_signal * crystal_gain_factor->GetBinContent(crystal + 1);
            mean_toa /= toa_used;
        }
        if (is_tot_event) {
            tot_events++;
            // continue;
        }

        float x_cog, y_cog;
        bool keep = calculate_cog(cog_distribution, signals);
        keep &= (!is_tot_event);
        if (!keep) {
            continue;
        }


        // Populate energies
        central_signal = signals[central_crystal_index];
        for (int i = 0; i < 8; i++) {
            int idx = center_nine_indexes[i];
            central_nine_signal += signals[idx];
        }
        for (int i = 0; i < sipms_to_use; i++) {
            int idx = remaining_indexes[i];
            total_signal += signals[idx];
        }
        if (signals[12] / total_signal < energy_fraction_cut) { // Skip events where the majority of the energy share is not in the central crystal;
            continue;
        }

        central_nine_signal += central_signal;
        total_signal += central_nine_signal;
        central_crystal_energy->Fill(central_signal);
        central_nine_energy->Fill(central_nine_signal);
        total_energy->Fill(total_signal);

        // Energy share fill
        for (int i = 0; i < 25; i++) {
            crystal_energy[i]->Fill(signals[i]);
            crystal_energy_shares[i]->Fill(signals[i] / total_signal);
        }

    }

    // fit_peak(central_crystal_energy);
    // fit_peak(central_nine_energy);
    // fit_peak(total_energy);


    int crystal_mapping[25] = {
        4, 9, 14, 19, 24,
        3, 8, 13, 18, 23,
        2, 7, 12, 17, 22,
        1, 6, 11, 16, 21,
        0, 5, 10, 15, 20
    };

    TCanvas* canvas = new TCanvas("gain_matching", "", 800, 600);
    canvas->SetRightMargin(0.05);
    central_crystal_energy->Draw("HIST e");
    // auto fit = central_crystal_energy->GetFunction("final_fit");
    // fit->Draw("same");
    // draw_text(fit, run_number, beam_energy);
    canvas->SaveAs("output/adc_calibration.pdf(");

    central_nine_energy->Draw("HIST e");
    // fit = central_nine_energy->GetFunction("final_fit");
    // fit->Draw("same");
    // draw_text(fit, run_number, beam_energy);
    canvas->SaveAs("output/adc_calibration.pdf");

    total_energy->Draw("HIST e");
    // fit = total_energy->GetFunction("final_fit");
    // fit->Draw("same");
    // draw_text(fit, run_number, beam_energy);
    canvas->SaveAs("output/adc_calibration.pdf");


    canvas->SetRightMargin(0.1);
    cog_distribution->Draw("COLZ");
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

    canvas->SaveAs("output/adc_calibration.pdf");

    canvas->Clear();
    gPad->SetLogz(0);
    pedestals->Draw("COLZ");
    pedestals->GetYaxis()->SetRangeUser(0, 200);
    canvas->SaveAs("output/adc_calibration.pdf");
    canvas->Clear();
    canvas->Divide(1, 1);
    TH1F* pedestals_width = new TH1F("pedestals_width", "Pedestal Width vs Channel;Channel;Width (ADC)", 72*8, 0, 72*8);
    for (int channel = 1; channel <= 72*8; channel++) {
        TH1D *proj = pedestals->ProjectionY("_py", channel, channel);
        if (proj && proj->GetEntries() > 0) {
            std::cout << "Channel " << channel << " pedestal width: " << proj->GetStdDev() << std::endl; 
            pedestals_width->SetBinContent(channel, proj->GetStdDev());
        }
    }
    pedestals_width->Draw("P");
    pedestals_width->GetYaxis()->SetRangeUser(0, 5);
    canvas->SaveAs("output/adc_calibration.pdf");


    canvas->Clear();
    canvas->Divide(5, 5);
    for (int i = 0; i < 25; i++) {
        canvas->cd(i + 1);
        crystal_energy[crystal_mapping[i]]->Draw("HIST e");
        // gPad->SetLogx();
    }

    canvas->SaveAs("output/adc_calibration.pdf");
    

    std::vector<float> gev_calib;
    canvas->Clear();
    canvas->Divide(5, 5);
    float mean_calib = 0;
    for (int i = 0; i < 25; i++) {
        canvas->cd(i + 1);
        crystal_energy_shares[crystal_mapping[i]]->Draw("HIST e");
        crystal_energy_shares[crystal_mapping[i]]->GetXaxis()->SetRangeUser(0.001, 1);
        gPad->SetLogx();
        float signal_for_1gev = crystal_energy[crystal_mapping[i]]->GetMean() / crystal_energy_shares[crystal_mapping[i]]->GetMean();
        float x_coord = 0.4;
        if (i == 12) {
            x_coord = 0.15;
        }
        TLatex *text = new TLatex();
        text->SetNDC();
        text->SetTextSize(0.04);
        text->DrawLatex(x_coord, 0.8, Form("Signal for 1 GeV: %.1f", signal_for_1gev));
        if (crystal_mapping[i] == 9) {
            continue;
        }
        gev_calib.push_back(signal_for_1gev);
        mean_calib += signal_for_1gev;
    }
    canvas->SaveAs("output/adc_calibration.pdf");
    mean_calib /= 24;   // Since we are excluding crystal 9
    std::cout << "1 GeV signal calibration: " << mean_calib << std::endl;

    TFile* calib_file = TFile::Open("output/adc_to_gev_calibration.root", "RECREATE");
    TParameter<float>* mean_calib_param = new TParameter<float>("mean_adc_to_gev_calibration", mean_calib);
    mean_calib_param->Write();
    calib_file->Close();


    for (int i = 0; i < 25; i++) {
        canvas->Clear();
        canvas->Divide(4, 4);
        for (int j = 0; j < sipms_to_use; j++) {
            canvas->cd(j + 1);
            sipm_energy[i][j]->Draw("HIST e");
            // canvas->GetXaxis()->SetMinimum(1);
            // gPad->SetLogx();
        }
        canvas->SaveAs("output/adc_calibration.pdf");
        canvas->Clear();
        canvas->Divide(4, 4);
        for (int j = 0; j < sipms_to_use; j++) {
            canvas->cd(j + 1);
            sipm_waveform[i][j]->Draw("COLZ");
            gPad->SetLogz();
        }
        canvas->SaveAs("output/adc_calibration.pdf");
        canvas->Clear();
        crystal_common_mode[i]->Draw("COLZ");
        gPad->SetLogz();
        canvas->SaveAs("output/adc_calibration.pdf");
    }
    gPad->SetLogz(0);
    
    canvas->Clear();
    toa_distribution->Draw("HIST");
    canvas->SaveAs("output/adc_calibration.pdf");

    canvas->Clear();
    toa_sample->Draw("HIST");
    canvas->SaveAs("output/adc_calibration.pdf");

    canvas->Clear();
    canvas->Divide(4, 4);
    auto text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    for (int sipm = 0; sipm < sipms_to_use; sipm++) {
        canvas->cd(sipm + 1);
        // auto fit = new TF1(Form("toa_fit_sipm_%d", sipm), "pol1", 0, 1024);
        E_vs_toa[sipm]->Draw("COLZ");
        // E_vs_toa[sipm]->Fit(fit, "RQ");
        // fit->SetLineColor(kRed);
        // text->DrawLatex(0.15, 0.85, Form("Slope: %.3f#pm%0.3f", fit->GetParameter(1), fit->GetParError(1)));
        E_vs_toa[sipm]->ProfileX()->Draw("same");

        gPad->SetLogz();
    }
    canvas->SaveAs("output/adc_calibration.pdf");
    gPad->SetLogz(0);

    auto canvas2 = new TCanvas("toa_correlations", "", 8000, 8000);
    canvas2->Divide(16, 16, 0.0005, 0.0005);
    for (int i = 0; i < sipms_to_use; i++) {
        for (int j = 0; j < sipms_to_use; j++) {
            if (j <= i) {
                canvas2->cd(i * 16 + j + 1);
                toa_correlations[i][j]->Draw("COLZ");
                // gPad->SetLogz();
            }
        }
    }
    canvas2->SaveAs("output/adc_correlation.png");

    canvas->Clear();
    max_sample_index->Draw("HIST");
    canvas->SaveAs("output/adc_calibration.pdf");

    canvas->Clear();
    second_max_sample_index->Draw("HIST");
    canvas->SaveAs("output/adc_calibration.pdf");

    canvas->Clear();
    third_max_sample_index->Draw("HIST");
    canvas->SaveAs("output/adc_calibration.pdf");


    canvas->Clear();
    canvas->SaveAs("output/adc_calibration.pdf)");

    float mean_x = cog_distribution->GetMean(1);
    float mean_y = cog_distribution->GetMean(2);
    float sigma_x = cog_distribution->GetStdDev(1);
    float sigma_y = cog_distribution->GetStdDev(2);

    std::cout << "Mean X: " << mean_x << " +/- " << sigma_x << std::endl;
    std::cout << "Mean Y: " << mean_y << " +/- " << sigma_y << std::endl;

    std::cout << "Total TOT events: " << tot_events << " out of " << nentries << std::endl;

}

