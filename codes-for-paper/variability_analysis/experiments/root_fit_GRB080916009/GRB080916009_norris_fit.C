// GRB080916009_norris_fit.C
//
// ROOT macro: read GRB080916009_lightcurve_10-400keV.csv (time_s, net_cts_per_s, net_err_cts_per_s --
// NaI detectors n3+n4 summed, background-subtracted, 10-400 keV) and fit a sum of 8 Norris (2005)
// pulses, same functional form and parameterization as this folder's Python pipeline
// (norris_fit.py::norris_pulse(), fitter.py's established 8-pulse decomposition: norris1/2/2-1/2-2/3/4/5/6).
//
// t_v itself (Bukhari et al. 2022 eq. 10) is NOT computed here -- this macro only reproduces the fit;
// t_v = (tau2/2)*sqrt((ln2+2*sqrt(tau1/tau2))^2 - 4*tau1/tau2) per pulse, same as the Python side.
//
// Starting parameters are the Python fit's converged values, rescaled from peak-normalized amplitude to
// this CSV's own physical counts/s scale (peak = 3241.92 cts/s) -- t_s/tau1/tau2 carry over unchanged,
// only A is rescaled. Refine further as needed; these are starting guesses, not a claim this CSV's own
// fit converges to the same numbers (different detector sum, different noise realization).
//
// Usage:
//   root -l 'GRB080916009_norris_fit.C("GRB080916009_lightcurve_10-400keV.csv")'
// or, from this folder with the default csv name:
//   root -l GRB080916009_norris_fit.C

#include <TAxis.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TMath.h>

#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

const int N_PULSES = 8;
const char *PULSE_NAMES[N_PULSES] = {
    "norris1", "norris2", "norris2-1", "norris2-2", "norris3", "norris4", "norris5", "norris6"
};

// Single Norris (2005) pulse, same functional form as norris_fit.py::norris_pulse():
//   A * exp(2*sqrt(tau1/tau2)) * exp(-tau1/(t-ts) - (t-ts)/tau2), t > ts, else 0.
double norris_pulse(double t, double A, double ts, double tau1, double tau2) {
    double dt = t - ts;
    if (dt <= 0) return 0.0;
    return A * TMath::Exp(2.0 * TMath::Sqrt(tau1 / tau2)) * TMath::Exp(-tau1 / dt - dt / tau2);
}

// Sum of N_PULSES Norris pulses. par layout: [A_i, ts_i, tau1_i, tau2_i] for i = 0..N_PULSES-1 (4*N_PULSES
// parameters total), same flat ordering as the Python fit's parameter vector.
double norris_sum(double *x, double *par) {
    double t = x[0];
    double total = 0.0;
    for (int i = 0; i < N_PULSES; i++) {
        total += norris_pulse(t, par[4 * i + 0], par[4 * i + 1], par[4 * i + 2], par[4 * i + 3]);
    }
    return total;
}

// use_errors = false: fit unweighted (TGraph, not TGraphErrors) -- ignores net_err_cts_per_s entirely.
// Added 2026-09-17 after the errors-weighted fit came back chi2/ndf ~ 0.03, evidence net_err_cts_per_s
// is inflated relative to the data's real point-to-point scatter (still unresolved which .dat column is
// responsible -- see the discussion this macro's first version prompted). Error bars are still read and
// available in yerr_vec below if you want to reintroduce them once that's sorted out.
void GRB080916009_norris_fit(const char *csv_path = "GRB080916009_lightcurve_10-400keV.csv", bool use_errors = true) {
    // --- Read CSV. Manual parse (not TTree::ReadFile) because the trailing "detectors_summed" column is
    // a string, which TTree::ReadFile's column auto-typing doesn't handle cleanly -- pulls only the three
    // numeric columns needed, by fixed index per this file's own header order:
    // time_s=0, net_cts_per_s=5, net_err_cts_per_s=6.
    std::ifstream file(csv_path);
    if (!file.is_open()) {
        std::cerr << "Could not open " << csv_path << std::endl;
        return;
    }
    std::string line;
    std::getline(file, line);  // header, discarded

    std::vector<double> t_vec, y_vec, yerr_vec;
    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::string field;
        std::vector<std::string> fields;
        while (std::getline(ss, field, ',')) fields.push_back(field);
        if (fields.size() < 7) continue;
        t_vec.push_back(std::stod(fields[0]));
        y_vec.push_back(std::stod(fields[5]));
        yerr_vec.push_back(std::stod(fields[6]));
    }
    file.close();

    int n = (int)t_vec.size();
    std::cout << "Read " << n << " points from " << csv_path << std::endl;
    if (n == 0) {
        std::cerr << "No data rows read -- check csv_path." << std::endl;
        return;
    }

    TGraph *gr = use_errors ? (TGraph *)new TGraphErrors(n, t_vec.data(), y_vec.data(), nullptr, yerr_vec.data())
                            : (TGraph *)new TGraph(n, t_vec.data(), y_vec.data());
    std::cout << "Fitting " << (use_errors ? "WEIGHTED by net_err_cts_per_s" : "UNWEIGHTED (net_err_cts_per_s ignored)")
              << std::endl;
    gr->SetTitle(
        "GRB080916C 10-400 keV NaI (n3+n4), background-subtracted;"
        "Time since trigger [s];Net count rate [cts/s]"
    );
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.4);
    gr->SetMarkerColor(kAzure + 2);
    gr->SetLineColor(kAzure + 2);

    // --- Total-fit function: sum of 8 Norris pulses, over the same [-1, 70]s window used throughout
    // this project's Python pipeline (fitter.py's START1/END1).
    double fit_min = -1.0, fit_max = 70.0;
    TF1 *f_total = new TF1("f_total", norris_sum, fit_min, fit_max, 4 * N_PULSES);

    // Starting parameters: the established Python fit's converged values (fitter.py's P0 comment),
    // amplitude rescaled from peak-normalized (Python fit peak=1) to this CSV's own peak (3241.92 cts/s).
    // t_s/tau1/tau2 are shape parameters, independent of the y-axis scale, so they carry over unchanged.
    double p0[4 * N_PULSES] = {
        1397.5, -0.7127, 2.2796, 0.4810,    // norris1: precursor
        2435.5, 0.2496, 0.9850, 6.7772,     // norris2: broad envelope/tail (the one overshadowing 3/4/5)
        81.0, 1.3187, 0.6471, 0.0334,       // norris2-1: residual-seeded, t_peak ~1.47s
        51.1, 2.2369, 0.0196, 0.2618,       // norris2-2: residual-seeded, t_peak ~2.31s
        837.4, 5.2910, 0.7529, 0.4451,      // norris3: the ~5.7-5.9s peak
        1147.9, 5.2102, 29.1785, 14.5375,   // norris4: broad TR3 pedestal
        489.9, 53.6618, 19.8173, 0.7360,    // norris5
        500.8, 61.5566, 0.4016, 2.8796,     // norris6
    };
    for (int i = 0; i < 4 * N_PULSES; i++) f_total->SetParameter(i, p0[i]);

    // Same box constraints as NorrisFitter.fit_boundaries() (norris_fit.py:49-50): A >= 0, t_s in
    // [fit_min, fit_max], tau1/tau2 >= 1e-4 (no upper bound there; capped here only to keep Minuit's
    // search space finite -- raise if a fit wants to push past it).
    for (int i = 0; i < N_PULSES; i++) {
        f_total->SetParLimits(4 * i + 0, 0.0, 1.0e6);          // A
        f_total->SetParLimits(4 * i + 1, fit_min, fit_max);    // t_s
        f_total->SetParLimits(4 * i + 2, 1.0e-4, 1.0e4);       // tau1
        f_total->SetParLimits(4 * i + 3, 1.0e-4, 1.0e4);       // tau2
    }
    for (int i = 0; i < N_PULSES; i++) {
        f_total->SetParName(4 * i + 0, Form("%s_A", PULSE_NAMES[i]));
        f_total->SetParName(4 * i + 1, Form("%s_ts", PULSE_NAMES[i]));
        f_total->SetParName(4 * i + 2, Form("%s_tau1", PULSE_NAMES[i]));
        f_total->SetParName(4 * i + 3, Form("%s_tau2", PULSE_NAMES[i]));
    }

    // --- Fit. R = use f_total's own range: S = return a TFitResultPtr; M = MINOS-improved errors.
    TFitResultPtr fit_result = gr->Fit(f_total, "RSM");
    fit_result->Print();

    std::cout << "\nchi2/ndf = " << f_total->GetChisquare() << " / " << f_total->GetNDF() << " = "
              << f_total->GetChisquare() / f_total->GetNDF() << std::endl;

    // --- Plot: data + total fit + each individual pulse, same visual convention as the Python side's
    // plot_fit(show_individuals=True).
    TCanvas *c = new TCanvas("c", "GRB080916C Norris fit", 1400, 800);
    gr->Draw("AP");
    f_total->SetLineColor(kBlack);
    f_total->SetLineWidth(2);
    f_total->Draw("SAME");

    int colors[N_PULSES] = {kOrange + 1, kGreen + 2, kSpring, kViolet, kRed, kMagenta, kCyan + 2, kPink + 7};
    TF1 *f_individual[N_PULSES];
    TLegend *leg = new TLegend(0.65, 0.55, 0.98, 0.92);
    leg->AddEntry(gr, use_errors ? "10-400 keV NaI (n3+n4), bkg-subtracted" : "10-400 keV NaI (n3+n4), bkg-subtracted (unweighted fit)", "lep");
    leg->AddEntry(f_total, "Total fit (8 Norris pulses)", "l");
    for (int i = 0; i < N_PULSES; i++) {
        f_individual[i] = new TF1(
            Form("f_%s", PULSE_NAMES[i]),
            [i](double *x, double *par) { return norris_pulse(x[0], par[0], par[1], par[2], par[3]); },
            fit_min, fit_max, 4
        );
        for (int j = 0; j < 4; j++) f_individual[i]->SetParameter(j, f_total->GetParameter(4 * i + j));
        f_individual[i]->SetLineColor(colors[i]);
        f_individual[i]->SetLineStyle(2);
        f_individual[i]->Draw("SAME");
        leg->AddEntry(f_individual[i], PULSE_NAMES[i], "l");
    }
    leg->Draw();
    const char *suffix = use_errors ? "weighted" : "unweighted";
    c->SaveAs(Form("GRB080916009_norris_fit_ROOT_%s.png", suffix));
    c->SaveAs(Form("GRB080916009_norris_fit_ROOT_%s.pdf", suffix));
    std::cout << "\nsaved GRB080916009_norris_fit_ROOT_" << suffix << ".png/.pdf" << std::endl;
}
