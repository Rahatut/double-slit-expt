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

    WaveFunction(double slit_distance, double screen_distance) 
        : d(slit_distance), distance_to_screen(screen_distance) {
        generateUnmeasuredWavefunction(); 
    }

    void generateUnmeasuredWavefunction() {
        values = linspace(-10.0, 10.0, 1000); // Create a range of angles

        probs.clear();
        for (double angle : values) {
            double path_difference = d * sin(angle);
            double intensity = pow(cos(M_PI * path_difference / distance_to_screen), 2);
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

    double evaluate(double angle) {
        // Calculate intensity based on the angle using the path difference
        double path_difference = d * sin(angle);
        return pow(cos(M_PI * path_difference / distance_to_screen), 2);
    }

    double evaluateNormalized(double angle) {
        return evaluate(angle) / norm; // Normalize using the precomputed norm
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

    DoubleSlitExperiment(double slit_distance, double screen_distance, int width, int height)
        : slit_dist(slit_distance), distance_to_screen(screen_distance), 
          screen_width(width), screen_height(height) {
        wavefunction = new WaveFunction(slit_dist, distance_to_screen);
    }

    ~DoubleSlitExperiment() {
        delete wavefunction;
    }

    void fireElectron() {
        double angle = wavefunction->measure();
        double detected_x = distance_to_screen * tan(angle);
        
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
            int screen_x = screen_width / 2 + detections_x[i] * (screen_width / 40); // Adjust scaling factor
            int screen_y = screen_height / 2 - detections_y[i]; // Invert Y to match screen coordinates
            
            Color color = (i % 2 == 0) ? BLUE : GREEN; 
            DrawCircle(screen_x, screen_y, 2, color); // Draw detected electrons
        }
    }

    void clearScreen() {
        detections_x.clear();
        detections_y.clear();
        delete wavefunction;
        wavefunction = new WaveFunction(slit_dist, distance_to_screen); // Reinitialize wavefunction
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
        double bin_height = (static_cast<double>(histogram[i]) / max_count) * (screen_height/1.5);
        // Draw the rectangle representing the bin
        DrawRectangle(i * bin_width, screen_height - bin_height, bin_width, bin_height, RED);
    }
}

};

int main() {
    const int screenWidth = 1600;
    const int screenHeight = 800;
    InitWindow(screenWidth, screenHeight, "Double Slit Experiment");

    DoubleSlitExperiment experiment(1.0, 10.0, screenWidth, screenHeight);
    experiment.resetExperiment(false); // Start with an unmeasured wavefunction
    experiment.fireElectronBeam(5000);  // Fire a beam of electrons

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);
        experiment.showScreen();  // Display the electron detections
        EndDrawing();
    }

    CloseWindow();

    // Create a separate window for the histogram
    InitWindow(screenWidth, screenHeight, "Histogram");

    while (!WindowShouldClose()) {
        BeginDrawing();
        ClearBackground(BLACK);
        experiment.drawHistogram(20); // Draw the histogram with 20 bins
        EndDrawing();
    }

    CloseWindow();
    return 0;
}
