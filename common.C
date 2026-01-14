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

float sigma_cut = 2;

float center_x = 1.95396;
float sigma_x = 0.184758 * sigma_cut;
float center_y = 1.9688;
float sigma_y = 0.195137 * sigma_cut;

// float adc_calib = 26704.4;  // 32444.1 Signal_ADC = 1 GeV 
// float adc_calib = 57854.3;  // v3

int signal_method = 2;


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

float calculate_signal_v2(uint32_t *adc_values, float gain) {
    // Pedestal is the mean of the first three samples
    float pedestal = (adc_values[0] + adc_values[1] + adc_values[2]) / 3.0f;
    
    // Signal is the sum of the three greatest values
    float signal = 0.0f;
    float max1 = 0.0f, max2 = 0.0f, max3 = 0.0f;
    
    for (int i = 3; i < 20; ++i) {
        float sample = adc_values[i] - pedestal;
        if (sample > max1) {
            max3 = max2;
            max2 = max1;
            max1 = sample;
        } else if (sample > max2) {
            max3 = max2;
            max2 = sample;
        } else if (sample > max3) {
            max3 = sample;
        }
    }
    signal = max1 + max2 + max3;
    return signal * gain;
}

float calculate_signal_v3(uint32_t *adc_values, float gain) {
    // Pedestal is the mean of the first three samples
    float pedestal = (adc_values[0] + adc_values[1] + adc_values[2]) / 3.0f;

    // Signal is the sum of all samples above pedestal
    float signal = 0.0f;
    for (int i = 3; i < 20; ++i) {
        float sample = adc_values[i] - pedestal;
        if (sample > 0) {
            signal += sample;
        }
    }
    return signal * gain;
}

float calculate_signal(uint32_t *adc_values, float gain) {
    if (signal_method == 2) {
        return calculate_signal_v2(adc_values, gain);
    }
    else if (signal_method == 3) {
        return calculate_signal_v3(adc_values, gain);
    }
}

bool calculate_signal(uint32_t *adc_values, uint32_t *tot_values, float gain, float &signal) {
    if (signal_method == 2) {
        // Check if there is a ToT value
        uint32_t tot = 0;
        for (int i = 0; i < 20; ++i) {
            if (tot_values[i] > 0) {
                tot = tot_values[i];
                break;
            }
        }
        // If there is no ToT value, return the single SiPM ADC signal
        if (tot == 0) {
            signal = calculate_signal(adc_values, gain);
            return false;
        }

        // Otherwise, just return the first ToT value
        signal = tot * gain;
        return true;
    }
    else if (signal_method == 3) {
        // Check if there is a ToT value
        uint32_t tot = 0;
        for (int i = 0; i < 20; ++i) {
            if (tot_values[i] > 0) {
                tot = tot_values[i];
                break;
            }
        }
        // If there is no ToT value, return the single SiPM ADC signal
        if (tot == 0) {
            signal = calculate_signal_v3(adc_values, gain);
            return false;
        }
        signal = tot;
        return true;
    }
    else {
        // Default to ADC signal
        signal = calculate_signal(adc_values, gain);
        return false;
    }
}

bool is_tot(uint32_t *tot_values) {
    for (int i = 0; i < 20; ++i) {
        // std::cout << "TOT Value[" << i << "]: " << tot_values[i] << std::endl;
        if (tot_values[i] > 50) {
            return true;
        }
    }
    return false;
}

int32_t get_toa(uint32_t *toa_values) {
    for (int i = 0; i < 20; ++i) {
        if (toa_values[i] > 0) {
            return toa_values[i];
        }
    }
    return -1;
}

void fit_peak(TH1* hist) {
    float mean = hist->GetMean();
    float rms = hist->GetRMS();
    float fit_min = mean - 1.5 * rms;
    float fit_max = mean + 1.5 * rms;
    TF1* rough_fit = new TF1("rough_fit", "gaus", fit_min, fit_max);
    hist->Fit(rough_fit, "RQ");
    float peak = rough_fit->GetParameter(1);
    float sigma = rough_fit->GetParameter(2);
    rough_fit->SetLineColor(kBlue);
    rough_fit->SetLineStyle(2);

    TF1* second_fit = new TF1("second_fit", "gaus", peak - sigma, peak + sigma);
    hist->Fit(second_fit, "RQ");
    peak = second_fit->GetParameter(1);
    sigma = second_fit->GetParameter(2);

    // TF1* final_fit = new TF1("final_fit", "crystalball", peak - sigma, peak + sigma);
    // final_fit->SetParameters(second_fit->GetParameter(0), peak, sigma, 1.5, 2.0);
    TF1* final_fit = new TF1("final_fit", "gaus", peak - sigma, peak + sigma);
    hist->Fit(final_fit, "R");
}

void draw_text(TF1* fit, int run_number, float beam_energy) {
    TLatex *text = new TLatex();
    text->SetNDC();
    text->SetTextSize(0.04);
    text->SetTextFont(42);
    text->SetTextAlign(31);

    float peak = fit->GetParameter(1);;
    float sigma = fit->GetParameter(2);
    float resolution = (sigma / peak) * 100.0f;
    text->DrawLatex(0.93, 0.85, Form("%.01f GeV Electrons", beam_energy));
    text->DrawLatex(0.93, 0.80, Form("Run %d", run_number));
    text->DrawLatex(0.93, 0.75, Form("Peak: %.03f", peak));
    text->DrawLatex(0.93, 0.70, Form("Sigma: %.03f", sigma));
    text->DrawLatex(0.93, 0.65, Form("Resolution: %.2f%%", resolution));
    text->DrawLatex(0.93, 0.60, Form("Signal method %d", signal_method));

}

bool position_cut(float x, float y) {
    // Cut anything outside of the center
    if (std::abs(x - center_x) > sigma_x || std::abs(y - center_y) > sigma_y) {
        return false;
    }
    return true;
}

bool calculate_cog(TH2* distribution, float *values) {
    float total_signal = 0;
    float x_weighted_sum = 0;
    float y_weighted_sum = 0;
    float x_cog = 0;
    float y_cog = 0;
    float w = 4;

    for (int i = 0; i < 25; i++) {
        total_signal += values[i];
    }
    
    float total_weight = 0;
    for (int i = 0; i < 25; i++) {
        int x = i % 5;
        int y = i / 5;
        float signal = values[i];
        float weight = w + std::log(signal / total_signal); // Avoid log(0)
        if (weight < 0) {
            weight = 0;
        }
        total_weight += weight;
        x_weighted_sum += x * weight;
        y_weighted_sum += y * weight;

    }
    if (total_weight > 0) {
        x_cog = x_weighted_sum / total_weight;
        y_cog = y_weighted_sum / total_weight;
        distribution->Fill(x_cog, y_cog);
    }
    return position_cut(x_cog, y_cog);
}

void print_progress(int progress) {
    std::cout << " [";
    for (int i = 0; i < 25; i++) {
        if (i < progress) {
            std::cout << "*";
        } else {
            std::cout << " ";
        }
    }
    std::cout << "]\r" << std::flush;
}