#include "footprint.h"
#include "input_read.h"
#include <atomic>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <m3dc1_source.h>
#include <memory>
#include <superposition_source.h>
#include <thread>
#include <unistd.h>

int main(int argc, char *argv[]) {
    // read params from input file
    std::string pathsFile = (argc > 1) ? argv[1] : "fpgen_input.txt";

    std::cout << "\n-----------------------------------------------\n"
              << "FPGEN - Magnetic Footprint Generator\n"
              << "-----------------------------------------------\n"
              << std::endl;
    input_read input(pathsFile);
    bool       readStatus = input.readInputFile();
    if (!readStatus) {
        std::cerr << "Error reading input file." << std::endl;
        return 1;
    }

    std::cout << "--------------- I/O PARAMETERS ----------------\n\n";
    std::cout << "nsources        = " << input.nsources << "\n";
    for (int i = 0; i < input.nsources; ++i)
        std::cout << "  source_" << i << "      = " << input.components[i].path
                  << "  ts=" << input.components[i].timeslice
                  << "  phase=" << input.components[i].phase << " deg"
                  << "  amp=" << input.components[i].amplitude << "\n";
    std::cout << "first_wall_path = " << input.first_wall_path << "\n"
              << "output_path     = " << input.output_path     << "\n\n"
              << "--------------- MAPPING PARAMETERS ------------\n\n"
              << "timeslice    = " << (input.nsources == 1
                                       ? std::to_string(input.components[0].timeslice)
                                       : "per-component") << "\n"
              << "manifold     = " << input.manifold    << "\n"
              << "grid_R1      = " << input.grid_R1     << "\n"
              << "grid_Z1      = " << input.grid_Z1     << "\n"
              << "grid_R2      = " << input.grid_R2     << "\n"
              << "grid_Z2      = " << input.grid_Z2     << "\n"
              << "nRZ          = " << input.nRZ         << "\n"
              << "nPhi         = " << input.nPhi        << "\n\n"
              << "--------------- INTEGRATOR PARAMETERS ---------\n\n"
              << "num_threads  = " << input.num_threads << "\n"
              << "max_turns    = " << input.max_turns   << "\n"
              << "h_init       = " << input.h_init      << "\n"
              << "h_min        = " << input.h_min       << "\n"
              << "h_max        = " << input.h_max       << "\n"
              << "-----------------------------------------------" << std::endl;

    // create footprint object
    std::cout << "\nCreating footprint object...\n"
              << std::endl;
    footprint footprint(input.manifold, input.grid_R1, input.grid_Z1, input.grid_R2, input.grid_Z2, input.nRZ, input.nPhi, input.max_turns);

    // set omp parameters
    omp_set_num_threads(input.num_threads);

    if (input.num_threads > 1) {
        std::cout << "\nA source and tracer object must be created for each thread. \nYou shall see the same message printed a number of times."
                  << std::endl;
    } else {
        std::cout << "\nCreating source and tracer objects..."
                  << std::endl;
    }

    // Each thread needs its own FieldSource (Fusion-IO is not thread-safe).
    // unique_ptr keeps each source at a stable address for the maglit references.
    auto make_source = [&]() -> std::unique_ptr<FieldSource> {
        if (input.nsources == 1) {
            return std::make_unique<M3DC1Source>(
                input.components[0].path.c_str(), input.components[0].timeslice);
        }
        auto ss = std::make_unique<SuperpositionSource>();
        for (const auto &comp : input.components)
            ss->add_component(comp.path, comp.timeslice,
                              comp.phase * M_PI / 180.0, comp.amplitude);
        return ss;
    };

    std::vector<std::unique_ptr<FieldSource>> sources;
    std::vector<std::unique_ptr<maglit>>      tracers;
    sources.reserve(input.num_threads);
    tracers.reserve(input.num_threads);
    for (int i = 0; i < input.num_threads; ++i) {
        sources.push_back(make_source());
        tracers.push_back(std::make_unique<maglit>(*sources.back()));
        tracers.back()->set_monitor(input.first_wall_path);
        tracers.back()->configure(input.h_init, input.h_min, input.h_max);
    }

    std::vector<maglit*> tracer_ptrs;
    tracer_ptrs.reserve(input.num_threads);
    for (auto &t : tracers)
        tracer_ptrs.push_back(t.get());

    std::cout << "\nSource and tracer object(s) created successfully."
              << std::endl;

    std::cout << "\nRunning grid...\n"
              << std::endl;

    // Per-thread progress counters for the multi-bar display.
    int total    = input.nPhi * input.nRZ;
    int nthreads = input.num_threads;
    int expected = (total + nthreads - 1) / nthreads;

    bool use_fancy = isatty(STDOUT_FILENO) && nthreads > 1;

    std::vector<footprint::ThreadProgress> thread_prog(nthreads);
    std::atomic<bool> disp_done{false};
    std::thread disp_thread;

    if (use_fancy) {
        const int W = 40;
        // Print initial 0% bars — these occupy the N lines the display thread will overwrite.
        for (int t = 0; t < nthreads; ++t)
            std::cout << "Thread " << std::setw(2) << t << " [>"
                      << std::string(W - 1, ' ') << "] 0/" << expected << "\n";
        std::cout.flush();

        disp_thread = std::thread([&, W]() {
            while (!disp_done.load(std::memory_order_relaxed)) {
                std::cout << "\033[" << nthreads << "A";
                for (int t = 0; t < nthreads; ++t) {
                    int done = thread_prog[t].value.load(std::memory_order_relaxed);
                    float frac = std::min(1.0f, (float)done / expected);
                    int filled = (int)(W * frac);
                    std::cout << "\rThread " << std::setw(2) << t << " ["
                              << std::string(filled, '=')
                              << (filled < W ? ">" : "")
                              << std::string(std::max(0, W - filled - 1), ' ')
                              << "] " << done << "/" << expected << "\n";
                }
                std::cout.flush();
                std::this_thread::sleep_for(std::chrono::milliseconds(150));
            }
            // Final frame — all bars filled
            std::cout << "\033[" << nthreads << "A";
            for (int t = 0; t < nthreads; ++t) {
                int done = thread_prog[t].value.load(std::memory_order_relaxed);
                std::cout << "\rThread " << std::setw(2) << t << " ["
                          << std::string(W, '=') << "] " << done << " done   \n";
            }
            std::cout.flush();
        });
    }

    footprint.run(tracer_ptrs, use_fancy ? &thread_prog : nullptr);

    if (use_fancy) {
        disp_done.store(true, std::memory_order_relaxed);
        disp_thread.join();
    }

    // Append .dat if no recognised extension is present
    auto has_ext = [](const std::string &p, const std::string &ext) {
        return p.size() >= ext.size() && p.substr(p.size() - ext.size()) == ext;
    };
    if (!has_ext(input.output_path, ".dat") &&
        !has_ext(input.output_path, ".txt") &&
        !has_ext(input.output_path, ".csv")) {
        input.output_path += ".dat";
    }

    // Collision check — append _new rather than overwriting
    std::string out_path = input.output_path;
    if (std::ifstream(out_path).is_open()) {
        auto dot = out_path.rfind('.');
        std::string base = (dot != std::string::npos) ? out_path.substr(0, dot) : out_path;
        std::string ext  = (dot != std::string::npos) ? out_path.substr(dot)    : ".dat";
        out_path = base + "_new" + ext;
        std::cerr << "Warning: output file exists. Saving to " << out_path
                  << " to avoid overwriting your data." << std::endl;
    }

    std::cout << "\nSaving output file at " << out_path << std::endl;
    if (!footprint.save(out_path)) {
        std::cerr << "Failed to save output file." << std::endl;
        return 1;
    }

    // print success message
    std::cout << "\nOutput file saved successfully.\n";
    std::cout << "\nProgram finished successfully.\n"
              << std::endl;

    return 0;
}