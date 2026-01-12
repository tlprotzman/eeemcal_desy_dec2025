#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TH2D.h>
#include <TCanvas.h>

std::map<int, std::vector<int>> read_mapping(const std::string& filename) {
    std::map<int, std::vector<int>> mapping;
    std::ifstream file(filename);
    std::string line;
    
    // Skip header
    std::getline(file, line);
    
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        int fpga, asic, connector, crystal, sipm, channel;
        char comma;
        
        // Parse CSV: FPGA,ASIC,Connector,Crystal,SiPM,Channel
        if (iss >> fpga >> comma >> asic >> comma >> connector >> comma >> crystal >> comma >> sipm >> comma >> channel) {
            // Ensure the vector is large enough for all SiPMs (16 per crystal)
            if (mapping[crystal].size() <= sipm) {
                mapping[crystal].resize(sipm + 1, -1);
            }
            mapping[crystal][sipm] = channel;
        }
    }
    
    file.close();
    return mapping;
}

void draw_waveform(int run_number) {
    char file_path[256];
    sprintf(file_path, "/Users/tristan/dropbox/eeemcal_desy_dec_2025/prod_0/Run%03d.root", run_number);
    
    auto mapping = read_mapping("eeemcal_desy_dec2025_mapping.csv");
    
    TFile* root_file = TFile::Open(file_path);
    TTree* tree = (TTree*)root_file->Get("events");
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("adc", 1);
    
    // Create histograms for all crystals and SiPMs
    std::map<std::string, TH2D*> histograms;
    
    for (int crystal = 0; crystal < 25; ++crystal) {
        auto it = mapping.find(crystal);
        if (it == mapping.end() || it->second.size() != 16) {
            if (it == mapping.end()) {
                std::cout << "No mapping for crystal " << crystal << std::endl;
            } else {
                std::cout << "Expected 16 channels for crystal " << crystal 
                          << ", found " << it->second.size() << std::endl;
            }
            continue;
        }
        
        for (int i = 0; i < 16; ++i) {
            char hist_name[64];
            char hist_title[128];
            sprintf(hist_name, "crystal_%d_sipm_%d_waveform", crystal, i);
            sprintf(hist_title, "Crystal %d SiPM %d Waveform;Time Bin;ADC", crystal, i);
            
            TH2D* hist = new TH2D(hist_name, hist_title, 20, 0, 20, 1024, 0, 1024);
            histograms[std::string(hist_name)] = hist;
        }
    }
    
    // Fill histograms from tree
    Long64_t nEntries = tree->GetEntries();
    int adc[576][20];
    tree->SetBranchAddress("adc", &adc);
    
    for (Long64_t entry = 0; entry < nEntries; ++entry) {
        tree->GetEntry(entry);
        
        for (int crystal = 0; crystal < 25; ++crystal) {
            auto it = mapping.find(crystal);
            if (it == mapping.end() || it->second.size() != 16) {
                continue;
            }
            
            const std::vector<int>& channels = it->second;
            for (int i = 0; i < 16; ++i) {
                char hist_name[64];
                sprintf(hist_name, "crystal_%d_sipm_%d_waveform", crystal, i);
                
                int ch = channels[i];
                TH2D* hist = histograms[std::string(hist_name)];
                
                for (int t = 0; t < 20; ++t) {
                    hist->Fill(t, adc[ch][t]);
                }
            }
        }
    }
    
    // Draw and save histograms
    TCanvas* canvas = new TCanvas("waveforms", "Waveforms", 1200, 800);
    
    char output_file[256];
    sprintf(output_file, "output/run%03d_waveform.pdf", run_number);
    
    for (int crystal = 0; crystal < 25; ++crystal) {
        auto it = mapping.find(crystal);
        if (it == mapping.end() || it->second.size() != 16) {
            continue;
        }
        
        canvas->Clear();
        canvas->Divide(4, 4);
        
        for (int i = 0; i < 16; ++i) {
            canvas->cd(i + 1);
            char hist_name[64];
            sprintf(hist_name, "crystal_%d_sipm_%d_waveform", crystal, i);
            
            auto hist_it = histograms.find(std::string(hist_name));
            if (hist_it != histograms.end()) {
                gPad->SetLogz();
                hist_it->second->Draw("COLZ");
            }
        }
        
        if (crystal == 0) {
            canvas->SaveAs(Form("%s(", output_file));
        } else if (crystal == 24){
            canvas->SaveAs(Form("%s)", output_file));
        } else {
            canvas->SaveAs(output_file);
        }
    }
    
    root_file->Close();
}
