#include <map>
#include <vector>

#include <TFitResult.h>
#include <TFitResultPtr.h>
#include <TF1.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TError.h>

#include "common.C"


void gain_match() {
    gStyle->SetOptStat(0);
    gErrorIgnoreLevel = kWarning;

    std::vector<int> run_numbers = {385, 386, 387, 388, 389, 390, 391, 392, 393, 394, 410, 396, 397, 398, 399, 400, 401, 402, 403, 404, 405, 406, 407, 408, 409};
    std::vector<int> run_crystal = { 19,  23,  24,  18,  13,  14,   8,   3,   4,   9,   2,   7,  12,  17,  22,  21,  16,  11,   6,   1,   0,   5,  10,  15,  20};
    std::vector<TH1*> histograms(400); // 25 crystals * 16 SiPMs

    std::cout << "matching " << run_crystal.size() << " crystals" << std::endl;

    // Set up histograms
    auto mapping = read_mapping("eeemcal_desy_dec2025_mapping.csv");
    for (int crystal = 0; crystal < 25; crystal++) {
        for (int sipm = 0; sipm < 16; sipm++) {
            int channel = mapping[crystal][sipm];
            auto hist = new TH1F(Form("crystal_%d_sipm_%d", crystal, sipm),
                                 Form("Crystal %d SiPM %d Signal;Signal (ADC);Counts", crystal, sipm),
                                 200, 0, 5500);
            histograms[crystal * 16 + sipm] = hist;
        }
    }


    // Process each run
    for (int i = 0; i < run_numbers.size(); i++) {
        print_progress(i);
        int run_number = run_numbers[i];
        int crystal = run_crystal[i];
        TFile* root_file = TFile::Open(Form("/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run%03d.root", run_number));
        TTree* tree = (TTree*)root_file->Get("events");
        uint32_t adc[576][20];
        tree->SetBranchAddress("adc", &adc);
        tree->SetBranchStatus("*", 0);
        tree->SetBranchStatus("adc", 1);

        // Fill histograms from tree
        Long64_t nEntries = tree->GetEntries();
        for (Long64_t entry = 0; entry < nEntries; ++entry) {
            tree->GetEntry(entry);
            for (int sipm = 0; sipm < 16; sipm++) {
                int channel = mapping[crystal][sipm];
                float signal = calculate_signal(adc[channel], 1.0f);
                histograms[crystal * 16 + sipm]->Fill(signal);
            }
        }
    }

    std::vector<float> peak_locations(400);

    // Fit each histogram to find the peak
    for (int channel = 0; channel < histograms.size(); channel++) {
        TH1* hist = histograms[channel];
        float mean = hist->GetMean();
        float rms = hist->GetRMS();
        float fit_min = mean - 1.5 * rms;
        float fit_max = mean + 1.5 * rms;
        TF1* rough_fit = new TF1("rough_fit", "gaus", fit_min, fit_max);
        auto fit_result = hist->Fit(rough_fit, "RQS");
        if (fit_result->Status() != 0) {
            std::cout << Form("Channel %d: Initial fit failed, skipping...", channel) << std::endl;
            continue;
            peak_locations[channel] = 0;
        }
        float peak = rough_fit->GetParameter(1);
        float sigma = rough_fit->GetParameter(2);
        rough_fit->SetLineColor(kBlue);
        rough_fit->SetLineStyle(2);

        TF1* second_fit = new TF1("second_fit", "gaus", peak - sigma, peak + sigma);
        hist->Fit(second_fit, "RQ");
        peak = second_fit->GetParameter(1);
        sigma = second_fit->GetParameter(2);

        TF1* final_fit = new TF1("final_fit", "crystalball", peak - sigma, peak + sigma);
        final_fit->SetParameters(second_fit->GetParameter(0), peak, sigma, 1.5, 2.0);
        hist->Fit(final_fit, "RMQ");

        peak_locations[channel] = final_fit->GetParameter(1);
    }

    // Print the average peak location
    float total_peak = 0.0f;
    int count = 0;
    for (float peak : peak_locations) {
        if (peak > 0) {
            total_peak += peak;
            count++;
        }
    }
    float average_peak = total_peak / count;
    std::cout << "Average Peak Location: " << average_peak << std::endl;

    // Calculate the mean signal per crystal
    std::vector<float> crystal_mean;
    for (int crystal = 0; crystal < 25; crystal++) {
        float mean_peak = 0;
        int ch_used = 0;
        for (int sipm = 0; sipm < 16; sipm++) {
            if ((crystal == 21) && (sipm == 1 || sipm == 4 || sipm == 10 || sipm == 11 || sipm == 13 || sipm == 14 || sipm == 15)) {    // Dead channels
                continue;
            }
            if (crystal == 15 && (sipm == 10 || sipm == 0 || sipm == 4 || sipm == 13 || sipm == 11)) {  // Looks very weird...
                continue;
            }
            mean_peak += peak_locations[crystal * 16 + sipm];
            ch_used++;
        }
        mean_peak /= ch_used;
        crystal_mean.push_back(mean_peak);
        std::cout << Form("Crystal %d peak: %.0f", crystal, mean_peak) << std::endl;
    }

    // Calculate the gain factors
    std::vector<float> gain_factors(400);
    if (false) {
        float target = 1300;//average_peak;
        for (int i = 0; i < peak_locations.size(); i++) {
            if (peak_locations[i] > 0) {
                gain_factors[i] = target / peak_locations[i];
            } else {
                gain_factors[i] = 1.0f; // Default factor if no peak
            }
        }
    }
    if (true) {
        for (int crystal = 0; crystal < 25; crystal++) {
            for (int sipm = 0; sipm < 16; sipm++) {
                int channel = crystal * 16 + sipm;
                if (peak_locations[channel] > 0) {
                    gain_factors[channel] = crystal_mean[crystal] / peak_locations[channel];
                } else {
                    gain_factors[channel] = 1;
                }
            }
        }
    }

    // Calculate per crystal gain factors
    float mean_crystal_signal = 0;
    for (int crystal = 0; crystal < 25; crystal++) {
        if (crystal == 9) {
            continue; // Bad crystal
        }
        mean_crystal_signal += crystal_mean[crystal];
    }
    mean_crystal_signal /= 24;
    std::cout << Form("Mean across all crystals: %.0f", mean_crystal_signal);
    std::vector<float> crystal_gain;
    for (int crystal = 0; crystal < 25; crystal++) {
        if (crystal_mean[crystal] > 0) {
            crystal_gain.push_back(mean_crystal_signal / crystal_mean[crystal]);
        } else {
            crystal_gain.push_back(1);
        }
    }

    gain_factors[15 * 16 + 10] = 0; // Crystal 15 SiPM 10
    gain_factors[21 * 16 + 1] = 0; // Crystal 21 SiPM 1
    gain_factors[21 * 16 + 4] = 0; // Crystal 21 SiPM 4
    gain_factors[21 * 16 + 10] = 0; // Crystal 21 SiPM 10
    gain_factors[21 * 16 + 11] = 0; // Crystal 21 SiPM 11
    gain_factors[21 * 16 + 13] = 0; // Crystal 21 SiPM 13
    gain_factors[21 * 16 + 14] = 0; // Crystal 21 SiPM 14
    gain_factors[21 * 16 + 15] = 0; // Crystal 21 SiPM 15

    // Save all histograms to a single PDF
    const char* output_file = "output/gain_matching.pdf";
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    text->SetTextFont(42);
    TCanvas* canvas = new TCanvas("gain_matching", "", 800, 600);
    for (int crystal = 0; crystal < run_crystal.size(); crystal++) {
        canvas->Clear();
        canvas->Divide(4, 4);
        for (int sipm = 0; sipm < 16; sipm++) {
            canvas->cd(sipm + 1);
            auto hist = histograms[crystal * 16 + sipm];
            hist->Draw();
            text->DrawLatex(0.15, 0.85, Form("Channel %d", crystal * 16 + sipm));
            text->DrawLatex(0.15, 0.80, Form("Ch gain: %.2f", gain_factors[crystal * 16 + sipm]));
            text->DrawLatex(0.15, 0.75, Form("Crystal gain: %.2f", crystal_gain[crystal]));
        }
        if (crystal == 0) {
            canvas->SaveAs(Form("%s(", output_file));
        } else {
            canvas->SaveAs(output_file);
        }
    }

    canvas->Clear();


    // Make a histogram of gain factors
    TH1F* gain_hist = new TH1F("gain_factors", "Gain Factors;Channel;Gain Factor", 400, 0, 400);
    for (int channel = 0; channel < gain_factors.size(); channel++) {
        gain_hist->SetBinContent(channel + 1, gain_factors[channel]);
        gain_hist->SetBinError(channel + 1, 0.0001);
        std::cout << "Channel " << channel << ": Gain Factor = " << gain_factors[channel] << std::endl;
    }

    gain_hist->SetMinimum(0);
    gain_hist->SetMaximum(2);
    gain_hist->Draw("e");
    canvas->SaveAs(Form("%s", output_file));
    
    TH1F* crystal_gain_hist = new TH1F("crystal_factor", "Crystal Gain;Crystal;Gain Factor", 25, 0, 25);
    for (int crystal = 0; crystal < 25; crystal++) {
        crystal_gain_hist->SetBinContent(crystal + 1, crystal_gain[crystal]);
        crystal_gain_hist->SetBinError(crystal + 1, 0.0001);
        std::cout << "Crystal " << crystal << ": Gain Factor = " << crystal_gain[crystal] << std::endl;
    }
    crystal_gain_hist->SetMinimum(0);
    crystal_gain_hist->SetMaximum(2);
    crystal_gain_hist->Draw("e");
    canvas->SaveAs(Form("%s)", output_file));

    
    
    TFile* out_root = new TFile("output/gain_factors.root", "RECREATE");
    gain_hist->Write();
    crystal_gain_hist->Write();
    out_root->Close();
}