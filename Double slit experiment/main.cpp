#include "raylib.h"
#include <vector>
#include <random>
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

class WaveFunction {
public:
    double d; // Slit distance
    double distance_to_screen; // Distance to the screen
    vector<double> values; 
    vector<double> probs; 
    double norm;  

    WaveFunction(double slit_distance, double screen_distance, bool measure_slit) 
        : d(slit_distance), distance_to_screen(screen_distance) {
        if (!measure_slit) {
            generateUnmeasuredWavefunction(); 
        } else {
            generateMeasuredWavefunction();
        }
    }

    void generateUnmeasuredWavefunction() {
        values = linspace(-10.0, 10.0, 1000); // Create a range of angles
        probs.clear();
        
        for (double angle : values) {
            double intensity = pow(cos(M_PI * d * sin(angle) / distance_to_screen), 2);
            probs.push_back(intensity);
        }

        // Normalize probabilities
        norm = trapezoidalRule(values);
        for (double& prob : probs) {
            prob /= norm; // Normalize probabilities
        }
    }

    void generateMeasuredWavefunction() {
        values.clear(); // Clear previous values
        values.push_back(-d / 2);
        values.push_back(d / 2);
        norm = 1;
        probs = {0.5, 0.5}; 
    }

    double trapezoidalRule(const vector<double>& points) {
        double sum = 0.0;
        size_t n = points.size();
        
        // Calculate the area using trapezoidal rule
        for (size_t i = 0; i < n - 1; ++i) {
            sum += (points[i] + points[i + 1]) / 2.0; // Average of the two points
        }

        double step_size = (points.back() - points.front()) / (n - 1); // Step size
        return sum * step_size; // Multiply by the width of the intervals
    }

    vector<double> linspace(double start, double end, int num) {
        vector<double> result;
        double step = (end - start) / (num - 1);
        for (int i = 0; i < num; ++i) {
            result.push_back(start + i * step);
        }
        return result;
    }

    double measure() {
        random_device rd;
        mt19937 gen(rd());
        discrete_distribution<> dist(probs.begin(), probs.end());
        return values[dist(gen)]; 
    }
};

class DoubleSlitExperiment {
public:
    double slit_dist; // Distance between slits
    double distance_to_screen; // Distance to screen
    int screen_width; // Width of the screen
    int screen_height; // Height of the screen
    vector<double> detections_x; // Detected x positions
    vector<double> detections_y; // Detected y positions
    WaveFunction* wavefunction; // Pointer to the wavefunction

    DoubleSlitExperiment(double slit_distance, double screen_distance, int width, int height, bool measure_slit)
        : slit_dist(slit_distance), distance_to_screen(screen_distance), 
          screen_width(width), screen_height(height) {
        wavefunction = new WaveFunction(slit_dist, distance_to_screen, measure_slit);
    }

    ~DoubleSlitExperiment() {
        delete wavefunction;
    }

    void fireElectron() {
        double angle = wavefunction->measure();
        double detected_x = distance_to_screen * sin(angle); // Corrected calculation
        
        detections_x.push_back(detected_x);
        detections_y.push_back(GetRandomValue(-screen_height / 2, screen_height / 2)); // Random y positions
    }

    void fireElectronBeam(int num_electrons) {
        for (int i = 0; i < num_electrons; ++i) {
            fireElectron();
        }
    }

    void showScreen() {
        for (size_t i = 0; i < detections_x.size(); ++i) {
            int screen_x = screen_width / 2 + static_cast<int>(detections_x[i] * (screen_width / 50)); // Adjust scaling factor
            int screen_y = screen_height / 2 - static_cast<int>(detections_y[i]); // Invert Y to match screen coordinates
            
            Color color = (i % 2 == 0) ? BLUE : GREEN; 
            DrawCircle(screen_x, screen_y, 1, color); // Draw detected electrons
        }
    }

    void clearScreen() {
        detections_x.clear();
        detections_y.clear();
        delete wavefunction;
        wavefunction = new WaveFunction(slit_dist, distance_to_screen, false); // Reinitialize wavefunction
    }

    void resetExperiment(bool measured) {
        if (measured) {
            wavefunction->generateMeasuredWavefunction();
        } else {
            wavefunction->generateUnmeasuredWavefunction();
        }
        clearScreen(); // Clear previous detections
    }

vector<int> calculateHistogram(int num_bins) {
        vector<int> histogram(num_bins, 0);
        double min_x = -screen_width / 2;
        double max_x = screen_width / 2;
        double bin_width = (max_x - min_x) / num_bins;

        for (double x : detections_x) {
            int bin_index = (x - min_x) / bin_width;
            if (bin_index >= 0 && bin_index < num_bins) {
                histogram[bin_index]++;
            }
        }
        return histogram;
    }

    void drawHistogram(int num_bins) {
    if (detections_x.empty()) {
        // If there are no detections, do not attempt to draw the histogram
        return;
    }

    vector<int> histogram = calculateHistogram(num_bins);
    
    // Find the maximum count for scaling
    int max_count = *max_element(histogram.begin(), histogram.end());

    // Handle the case where max_count is 0 to avoid division by zero
    if (max_count == 0) {
        return; // If no counts, do not draw the histogram
    }

    double bin_width = (screen_width) / num_bins;

    for (int i = 0; i < num_bins; ++i) {
        // Calculate the height of each bin
        double bin_height = (static_cast<double>(histogram[i]) / max_count) * (screen_height/1.2);
        // Draw the rectangle representing the bin
        DrawRectangle(i * bin_width, screen_height - bin_height, bin_width, bin_height, RED);
    }
}

    vector<int> findPeaks(const vector<int>& histogram, double threshold) {
        vector<int> peaks;
        int n = histogram.size();

        for (int i = 1; i < n - 1; ++i) {
            // Check if the current bin is a peak
            if (histogram[i] > threshold && histogram[i] > histogram[i - 1] && histogram[i] > histogram[i + 1]) {
                peaks.push_back(i);
            }
        }

        return peaks;
    }

    double averageDistanceBetweenPeaks(const vector<int>& peaks, double bin_width) {
        if (peaks.size() < 2) return 0.0; // Not enough peaks to calculate distance

        double total_distance = 0.0;
        for (size_t i = 1; i < peaks.size(); ++i) {
            total_distance += (peaks[i] - peaks[i - 1]) * bin_width; // Calculate distance in terms of x-coordinates
        }

        return total_distance / (peaks.size() - 1); // Return average distance
    }

    void analyzePeaks() {
        vector<int> histogram = calculateHistogram(20); // Number of bins for the histogram
        int max_count = *max_element(histogram.begin(), histogram.end());
        double threshold = static_cast<double>(max_count) * 0.5; // Threshold to identify peaks (50% of max count)
        
        vector<int> peaks = findPeaks(histogram, threshold);
        double bin_width = static_cast<double>(screen_width) / 20; // Width of each bin

        double bright_spot_distance = averageDistanceBetweenPeaks(peaks, bin_width);
        
    }
};

int main() {
    const int screenWidth = 1000;
    const int screenHeight = 600;
    const float displayTime = 0.05f;
    
    // Initialize the main window for the double slit experiment
    InitWindow(screenWidth, screenHeight, "Double Slit Experiment");
    
    // Loop through different slit distances
    for (double screen_distance = 10.0; screen_distance <= 30.0; screen_distance += 2.0) {
        DoubleSlitExperiment experiment(20.0, screen_distance, screenWidth, screenHeight, false); // Start with an unmeasured wavefunction
        experiment.fireElectronBeam(50000);  // Fire a beam of electrons

        // Analyze peaks after firing the electron beam
        experiment.analyzePeaks();

        // Display the results for the specified duration
        float startTime = GetTime();
        while (GetTime() - startTime < displayTime) {
            BeginDrawing();
            ClearBackground(BLACK);
            experiment.showScreen();  // Display the electron detections
            EndDrawing();
        }

        experiment.clearScreen(); // Clear detections for the next iteration
    }
    // Create a separate window for the histogram
    //InitWindow(screenWidth, screenHeight, "Histogram");

    //while (!WindowShouldClose()) {
       // BeginDrawing();
       // ClearBackground(BLACK);
       // experiment.drawHistogram(20); // Draw the histogram with 10 bins
        //EndDrawing();
   // }

   // CloseWindow();
    return 0;
}
