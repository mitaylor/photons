#ifndef PAPER_H
#define PAPER_H

#include <array>
#include <cmath>
#include <functional>
#include <string>
#include <vector>

#include "TCanvas.h"
#include "TObject.h"
#include "TGraph.h"
#include "TH1.h"
#include "TLegend.h"

#include "pencil.h"

class paper {
  public:
    paper(std::string const& tag, pencil* p, int64_t cols, int64_t rows)
        : _tag(tag),
          _size(0),
          _cols(cols),
          _rows(rows),
          _pencil(p),
          canvas(nullptr) { }

    paper(std::string const& tag)
        : paper(tag, 0, 0, 0) { }

    paper(std::string const& tag, pencil* p)
        : paper(tag, p, 0, 0) { }

    paper(std::string const& tag, int64_t cols, int64_t rows)
        : paper(tag, 0, cols, rows) { }

    paper(paper const&) = delete;
    paper& operator=(paper const&) = delete;
    ~paper() = default;

    void add() { ++_size; }

    void stack(int64_t index, TObject* const object) {
        objects.push_back(object);
        indices.push_back(index);
    }

    template <typename... T>
    void stack(int64_t index, TObject* const object, T const&... adjectives) {
        stack(index, object); _pencil->describe(object, adjectives...); }

    void add(TObject* const object) { add(); stack(_size, object); }

    template <typename... T>
    void add(TObject* const object, T const&... adjectives) {
        add(object); _pencil->describe(object, adjectives...); }

    void stack(TObject* const object) { stack(_size, object); }

    template <typename... T>
    void stack(TObject* const object, T const&... adjectives) {
        stack(object); _pencil->describe(object, adjectives...); }

    void divide(int64_t cols, int64_t rows) { _cols = cols; _rows = rows; }

    void decorate(std::function<void()> d) { _d = d; }
    void format(std::function<void(TH1*)> f) { _f = f; }
    void format(std::function<void(TGraph*)> g) { _g = g; }
    void legend(std::function<std::array<float, 4>()> l) { _l = l; }
    void style(std::function<void(TLegend*)> s) { _s = s; }

    void describe(pencil* pencil) { _pencil = pencil; }

    void draw(char const* ext) {
        using namespace std::literals::string_literals;

        if (_cols * _rows == 0) { split(); }

        if (canvas == NULL) {
            canvas = new TCanvas(("paper_"s + _tag).data(), "",
                             400 * _cols, 400 * _rows);
            canvas->Divide(_cols, _rows);

            int64_t count = static_cast<int64_t>(objects.size());
            for (int64_t i = 0; i < count; ++i) {
                canvas->cd(indices[i]);
                apply(objects[i], _f);
                apply(objects[i], _g);
                objects[i]->Draw("same p e");
                apply(_d);
            }

            legends();
        }

        save(ext);
    }

  private:
    template <typename T>
    void apply(TObject* const obj, std::function<void(T*)> f) {
        if (f && obj->InheritsFrom(T::Class()))
            f(static_cast<T*>(obj));
    }

    template <typename T>
    void apply(std::function<T> f) { if (f) { f(); } }

    void split() {
        float rows = std::ceil(std::sqrt(_size));
        float cols = std::ceil(_size / rows);

        _cols = cols;
        _rows = rows;
    }

    void legends() {
        if (_pencil == nullptr) { return; }

        auto description = _pencil->description();

        for (int64_t i = 0; i < _size; ++i) {
            int64_t index = i + 1;

            std::vector<TObject*> associates;
            int64_t count = static_cast<int64_t>(objects.size());
            for (int64_t j = 0; j < count; ++j)
                if (indices[j] == index)
                    associates.push_back(objects[j]);

            auto point = _l ? _l() : std::array<float, 4>{
                0.5, 0.9, 0.87, 0.04 };
            point[3] = point[2] - associates.size() * point[3];

            canvas->cd(index);

            TLegend* l = new TLegend(point[0], point[3], point[1], point[2]);
            apply(l, _s);

            for (auto const& obj : associates) {
                auto desc = description.find(obj) != std::end(description) ?
                    description[obj] : std::string(obj->GetName());
                l->AddEntry(obj, desc.data(), "pl");
            }

            l->Draw();
        }
    }

    void save(char const* ext) {
        using namespace std::literals::string_literals;

        canvas->SaveAs((_tag + "."s + ext).data());
    }

    std::string const _tag;
    int64_t _size;
    int64_t _cols;
    int64_t _rows;

    std::vector<TObject*> objects;
    std::vector<int64_t> indices;

    std::function<void()> _d;
    std::function<void(TH1*)> _f;
    std::function<void(TGraph*)> _g;
    std::function<std::array<float, 4>()> _l;
    std::function<void(TLegend*)> _s;

    pencil* _pencil;

    TCanvas* canvas;
};

#endif /* PAPER_H */
