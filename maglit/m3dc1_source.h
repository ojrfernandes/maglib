#pragma once
// Last modified: 26.05.30

#include "field_source.h"
#include "fusion_io.h"

// FieldSource implementation backed by an M3D-C1 HDF5 file via Fusion-IO.
//
// The Fusion-IO search hint is allocated in the constructor and freed in the
// destructor, so callers do not need to manage hint lifetime manually.
// psin and psi fields are loaded lazily on first use.
class M3DC1Source : public FieldSource {
  public:
    M3DC1Source(const char *path, int timeslice);
    ~M3DC1Source() override;

    // Non-copyable: owns Fusion-IO file handles.
    M3DC1Source(const M3DC1Source &)            = delete;
    M3DC1Source &operator=(const M3DC1Source &) = delete;

    bool eval_B(double R, double phi, double Z, double B[3]) override;
    bool eval_psin(double R, double phi, double Z, double &psin) override;
    bool eval_psi(double R, double phi, double Z, double &psi) override;

    // Returns false if the source or magnetic field failed to open.
    bool is_valid() const { return src_ != nullptr && b_field_ != nullptr; }

  private:
    fio_source      *src_         = nullptr;
    fio_field       *b_field_     = nullptr;
    fio_field       *psin_field_  = nullptr;
    fio_field       *psi_field_   = nullptr;
    fio_option_list  opt_;
    fio_hint         hint_        = nullptr;
};
