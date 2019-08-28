#include "../include/lambdas.h"

#include "../git/config/include/configurer.h"

#include "../git/history/include/interval.h"
#include "../git/history/include/multival.h"
#include "../git/history/include/history.h"

#include "../git/paper-and-pencil/include/paper.h"
#include "../git/paper-and-pencil/include/pencil.h"

#include "../git/tricks-and-treats/include/trunk.h"
#include "../git/tricks-and-treats/include/zip.h"

#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TLatex.h"

#include <functional>
#include <string>
#include <vector>

using namespace std::literals::string_literals;
using namespace std::placeholders;

int obnubilate(char const* config, char const* output) {
    auto conf = new configurer(config);

    auto ref = conf->get<std::string>("ref");
    auto label = conf->get<std::string>("label");
    auto tag = conf->get<std::string>("tag");
    auto system = conf->get<std::string>("system");

    auto inputs = conf->get<std::vector<std::string>>("inputs");
    auto labels = conf->get<std::vector<std::string>>("labels");
    auto legends = conf->get<std::vector<std::string>>("legends");
    auto figures = conf->get<std::vector<std::string>>("figures");
    auto columns = conf->get<std::vector<int32_t>>("columns");
    auto ranges = conf->get<std::vector<float>>("ranges");
    auto groups = conf->get<std::vector<int32_t>>("groups");

    auto info = conf->get<std::string>("info");

    auto dpt = conf->get<std::vector<float>>("pt_diff");
    auto dhf = conf->get<std::vector<float>>("hf_diff");
    auto dcent = conf->get<std::vector<int32_t>>("cent_diff");

    /* manage memory manually */
    TH1::AddDirectory(false);
    TH1::SetDefaultSumw2();

    /* open input files */
    TFile* f = new TFile(ref.data(), "read");

    std::vector<TFile*> files(inputs.size(), nullptr);
    zip([&](auto& file, auto const& input) {
        file = new TFile(input.data(), "read");
    }, files, inputs);

    /* prepare plots */
    auto hb = new pencil();
    hb->category("type", "total", labels);

    zip([&](auto const& label, auto const& legend) {
        hb->alias(label, legend); }, labels, legends);

    auto cs = std::vector<paper*>(figures.size(), nullptr);

    /* lambdas */
    std::function<void(TH1*)> square_ = [](TH1* h) {
        for_contents([](std::array<double, 1> v) {
            return v[0] * v[0]; }, h); };
    std::function<void(TH1*)> sqrt_ = [](TH1* h) {
        for_contents([](std::array<double, 1> v) {
            return std::sqrt(v[0]); }, h); };

    auto shader = [&](TH1* h, float max) {
        default_formatter(h, 0., max);
        auto col = h->GetLineColor();
        h->SetFillColorAlpha(col, 0.32);
        h->SetLineWidth(1);
    };

    auto pt_info = [&](int64_t index, float pos) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%.0f < p_{T}^{#gamma} < %.0f",
            dpt[index - 1], dpt[index]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, pos, buffer);
    };

    auto hf_info = [&](int64_t index, float pos) {
        char buffer[128] = { '\0' };
        sprintf(buffer, "%i - %i%%", dcent[index], dcent[index - 1]);

        TLatex* l = new TLatex();
        l->SetTextFont(43);
        l->SetTextSize(12);
        l->DrawLatexNDC(0.135, pos, buffer);
    };

    auto info_text = [&](int64_t index, history* h) {
        auto indices = h->indices_for(index - 1);
        auto pt_x = indices[0];
        auto hf_x = indices[1];

        pt_info(pt_x + 1, 0.75);
        hf_info(hf_x + 1, 0.71);
    };

    /* prepare output */
    TFile* fout = new TFile(output, "recreate");

    /* calculate variations */
    zip([&](auto const& figure, auto cols, auto range, auto& c) {
        auto stub = "_"s + figure;

        c = new paper(tag + "_var"s + stub, hb);
        apply_style(c, system + " #sqrt{s_{NN}} = 5.02 TeV",
            std::bind(shader, _1, range));
        c->divide(cols, -1);

        auto base = new history(f, tag + "_"s + label + stub, "base");

        std::vector<history*> sets;

        std::vector<history*> batches(inputs.size(), nullptr);
        zip([&](auto& batch, auto file, auto const& label) {
            batch = new history(file, tag + "_"s + label + stub, "batch");
        }, batches, files, labels);

        auto total = new history(*base, "total");
        total->apply([](TH1* h) { h->Reset("MICES"); });

        for (auto const& batch : batches) {
            batch->add(*base, -1);
            batch->apply(square_);
        }

        zip([&](auto const& batch, auto group) {
            if (group == static_cast<int32_t>(sets.size())) {
                sets.push_back(new history(*batch, "set"));
            } else {
                sets[group]->apply([&](TH1* h, int64_t i) {
                    for_contents([](std::array<double, 2> v) -> float {
                        return std::max(v[0], v[1]);
                    }, h, (*batch)[i]); });
            }

            batch->apply(sqrt_);
        }, batches, groups);

        for (auto const& set : sets)
            total->add(*set, 1);

        total->apply(sqrt_);

        /* add plots */
        auto style = [&](TH1* h) { c->adjust(h, "hist", "f"); };

        total->apply([&](TH1* h) { c->add(h, "total"); style(h); });
        zip([&](auto& batch, auto const& label) {
            batch->apply([&](TH1* h, int64_t index) {
                c->stack(index + 1, h, label); style(h);
            });

            batch->save(tag);
        }, batches, labels);

        /* add info text */
        if (info == "pt") { c->accessory(std::bind(pt_info, _1, 0.75)); }
        if (info == "hf") { c->accessory(std::bind(hf_info, _1, 0.75)); }
        if (info == "pthf") { c->accessory(std::bind(info_text, _1, base)); }

        /* save histograms */
        for (auto const& set : sets)
            set->save(tag);

        base->save(tag);
        total->save(tag);
    }, figures, columns, ranges, cs);

    /* draw plots */
    hb->sketch();
    for (auto c : cs)
        c->draw("pdf");

    fout->Close();

    return 0;
}

int main(int argc, char* argv[]) {
    if (argc == 3)
        return obnubilate(argv[1], argv[2]);

    return 0;
}