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

#include "common.C"

float energy_fraction_cut = 0.3;
long n_events = 1000000;
// n_events = 10000;

int run_number = 411;
float beam_energy = 1;


void adc_calibration() {
    gStyle->SetOptStat(0);
    float energy = 1;

    int central_crystal_index = 12;
    int center_nine_indexes[8] = {7, 8, 9, 11, 13, 17, 18, 19};
    int remaining_indexes[16] = {0, 1, 2, 3, 4, 5, 6, 10, 11, 15, 16, 20, 21, 22, 23, 24};


    TH1* central_crystal_energy = new TH1F("central_crystal_energy", "Central Crystal Energy;Energy (ADC);Events", 500, 0, 75000);
    TH1* central_nine_energy    = new TH1F("central_nine_energy", "Central 3x3 Energy;Energy (ADC);Events", 500, 0, 75000);
    TH1* total_energy           = new TH1F("total_energy", "Total Energy;Energy (ADC);Events", 500, 0, 75000);
    TH2* cog_distribution       = new TH2F("cog_distribution", "Center of Gravity Distribution;X (# Crystals);Y (# Crystals)", 100, -0.5, 4.5, 100, -0.5, 4.5);
    TH1 *toa_distribution        = new TH1F("toa_distribution", "ToA Distribution;ToA;Events", 1024, 0, 1024);
    
    std::vector<TH2*> E_vs_toa;
    for (int sipm = 0; sipm < 16; sipm++) {
        E_vs_toa.push_back(new TH2F(Form("E_vs_toa_sipm_%d", sipm), Form("SiPM %d Energy vs ToA;ToA;Energy (ADC)", sipm), 1024, 0, 1024, 500, 0, 3000));
    }   
    
    std::vector<std::vector<TH2*>> toa_correlations;
    for (int i = 0; i < 16; i++) {
        toa_correlations.push_back(std::vector<TH2*>());
        for (int j = 0; j < 16; j++) {
            toa_correlations[i].push_back(new TH2F(Form("toa_correlation_sipm_%02d_sipm_%02d", i, j),
                                                  Form("SiPM %02d vs SiPM %02d ToA;SiPM %02d ToA;SiPM %02d ToA", i, j, i, j),
                                                  1024, 0, 1024, 1024, 0, 1024));
        }
    }

    std::vector<std::vector<TH1*>> sipm_energy;
    std::vector<TH1*> crystal_energy;
    std::vector<TH1*> crystal_energy_shares;
    for (int i = 0; i < 25; i++) {
        for (int j = 0; j < 16; j++) {
            sipm_energy.push_back(std::vector<TH1*>());
            sipm_energy[i].push_back(new TH1F(Form("crystal_%02d_sipm_%02d_energy", i, j), Form("Crystal %02d SiPM %02d Energy;Energy (ADC);Events", i, j), 500, 0, 4000));
        }   
        crystal_energy.push_back(new TH1F(Form("crystal_%02d_energy", i), Form("Crystal %02d Energy;Energy (ADC);Events", i), 500, 0, 40000));
        crystal_energy_shares.push_back(new TH1F(Form("crystal_%02d_energy_share", i), Form("Crystal %02d Energy Share;Energy Share;Events", i), 10000, 0, 1));
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
    char file_path[256];
    sprintf(file_path, "/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run%03d.root", run_number);
    TFile* root_file = TFile::Open(file_path);
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

        float signals[25];
        for (int crystal = 0; crystal < 25; crystal++) {
            float mean_toa = 0;
            uint32_t toa_used = 0;
            if (crystal == 9) {
                continue;
            }
            float crystal_signal = 0.0f;
            for (int sipm = 0; sipm < 16; sipm++) {
                int channel = mapping[crystal][sipm];
                float gain = gain_factor->GetBinContent(crystal * 16 + sipm + 1);
                float channel_signal = calculate_signal(adc[channel], gain);
                sipm_energy[crystal][sipm]->Fill(channel_signal);
                pedestals->Fill(channel, (adc[channel][0] + adc[channel][1] + adc[channel][2]) / 3.0f);
                crystal_signal += channel_signal;
                uint32_t this_toa = get_toa(toa[channel]);
                if (this_toa >= 0) {
                    if (crystal == 12) {
                        toa_distribution->Fill(this_toa);
                        E_vs_toa[sipm]->Fill(this_toa, channel_signal);
                        for (int other_sipm = 0; other_sipm < 16; other_sipm++) {
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
        for (int i = 0; i < 16; i++) {
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

    fit_peak(central_crystal_energy);
    fit_peak(central_nine_energy);
    fit_peak(total_energy);


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
    auto fit = central_crystal_energy->GetFunction("final_fit");
    fit->Draw("same");
    draw_text(fit, run_number, beam_energy);
    canvas->SaveAs("output/adc_calibration.pdf(");

    central_nine_energy->Draw("HIST e");
    fit = central_nine_energy->GetFunction("final_fit");
    fit->Draw("same");
    draw_text(fit, run_number, beam_energy);
    canvas->SaveAs("output/adc_calibration.pdf");

    total_energy->Draw("HIST e");
    fit = total_energy->GetFunction("final_fit");
    fit->Draw("same");
    draw_text(fit, run_number, beam_energy);
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
    pedestals->Draw("COLZ");
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

    for (int i = 0; i < 25; i++) {
        canvas->Clear();
        canvas->Divide(4, 4);
        for (int j = 0; j < 16; j++) {
            canvas->cd(j + 1);
            sipm_energy[i][j]->Draw("HIST e");
            // canvas->GetXaxis()->SetMinimum(1);
            // gPad->SetLogx();
        }
        canvas->SaveAs("output/adc_calibration.pdf");
    }
    
    canvas->Clear();
    toa_distribution->Draw("HIST");
    canvas->SaveAs("output/adc_calibration.pdf");

    canvas->Clear();
    canvas->Divide(4, 4);
    auto text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    for (int sipm = 0; sipm < 16; sipm++) {
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
    for (int i = 0; i < 16; i++) {
        for (int j = 0; j < 16; j++) {
            if (j <= i) {
                canvas2->cd(i * 16 + j + 1);
                toa_correlations[i][j]->Draw("COLZ");
                // gPad->SetLogz();
            }
        }
    }
    canvas2->SaveAs("output/adc_correlation.png");

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

