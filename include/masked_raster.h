#pragma once
#include "pcrtypes.h"
#include "raster.h"


template<
    typename T>
class MaskedRaster:
    public Raster<T>
{

public:

                   MaskedRaster        ();

                   MaskedRaster        (size_t nr_rows,
                                        size_t nr_cols);

                   MaskedRaster        (MaskedRaster const& other)=default;

                   MaskedRaster        (MaskedRaster&& other)=default;

  virtual          ~MaskedRaster       ()=default;

  MaskedRaster&    operator=           (MaskedRaster const& other)=default;

  MaskedRaster&    operator=           (MaskedRaster&& other)=default;

  bool             is_mv               (size_t index) const;

  bool             is_mv               (size_t row,
                                        size_t col) const;

  void             set_mv              (size_t index);

  void             set_mv              (size_t row,
                                        size_t col);

private:

};


template<
    typename T>
MaskedRaster<T>::MaskedRaster()

    : Raster<T>()

{
}


template<
    typename T>
MaskedRaster<T>::MaskedRaster(
    size_t nr_rows,
    size_t nr_cols)

    : Raster<T>(nr_rows, nr_cols)

{
}


template<
    typename T>
bool MaskedRaster<T>::is_mv(
    size_t index) const
{
    return pcr::isMV<T>(this->cell(index));
}


template<
    typename T>
bool MaskedRaster<T>::is_mv(
    size_t row,
    size_t col) const
{
    return pcr::isMV<T>(this->cell(row, col));
}


template<
    typename T>
void MaskedRaster<T>::set_mv(
    size_t index)
{
    pcr::setMV<T>(this->cell(index));
}


template<
    typename T>
void MaskedRaster<T>::set_mv(
    size_t row,
    size_t col)
{
    pcr::setMV<T>(this->cell(row, col));
}
