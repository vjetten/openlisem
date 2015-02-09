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

                   MaskedRaster        (std::initializer_list<
                                            std::initializer_list<T>> const&
                                                values);

                   MaskedRaster        (size_t nr_rows,
                                        size_t nr_cols,
                                        double north=0.0,
                                        double west=0.0,
                                        double cell_size=1.0);

                   MaskedRaster        (MaskedRaster const& other)=default;

                   MaskedRaster        (MaskedRaster&& other)=default;

    virtual        ~MaskedRaster       ()=default;

    MaskedRaster&  operator=           (MaskedRaster const& other)=default;

    MaskedRaster&  operator=           (MaskedRaster&& other)=default;

    bool           is_mv               (size_t index) const;

    bool           is_mv               (size_t row,
                                        size_t col) const;

    void           set_mv              (size_t index);

    void           set_mv              (size_t row,
                                        size_t col);

    void           set_all_mv          ();

    void           replace_with_mv     (T const& value);

private:

};


template<
    typename T>
inline MaskedRaster<T>::MaskedRaster()

    : Raster<T>()

{
}


template<
    typename T>
inline MaskedRaster<T>::MaskedRaster(
    std::initializer_list<std::initializer_list<T>> const& values)

    : Raster<T>(values)

{
}


template<
    typename T>
inline MaskedRaster<T>::MaskedRaster(
    size_t nr_rows,
    size_t nr_cols,
    double north,
    double west,
    double cell_size)

    : Raster<T>(nr_rows, nr_cols, north, west, cell_size)

{
}


template<
    typename T>
inline bool MaskedRaster<T>::is_mv(
    size_t index) const
{
    return pcr::isMV<T>(this->cell(index));
}


template<
    typename T>
inline bool MaskedRaster<T>::is_mv(
    size_t row,
    size_t col) const
{
    return pcr::isMV<T>(this->cell(row, col));
}


template<
    typename T>
inline void MaskedRaster<T>::set_mv(
    size_t index)
{
    pcr::setMV<T>(this->cell(index));
}


template<
    typename T>
inline void MaskedRaster<T>::set_mv(
    size_t row,
    size_t col)
{
    pcr::setMV<T>(this->cell(row, col));
}


template<
    typename T>
inline void MaskedRaster<T>::set_all_mv()
{
    pcr::setMV<T>(&this->cell(0), this->nr_cells());
}


template<
    typename T>
inline void MaskedRaster<T>::replace_with_mv(
    T const& value)
{
    T* it = &this->cell(0);

    for(size_t i = 0; i < this->nr_cells(); ++i, ++it) {
        if(*it == value) {
            pcr::setMV(*it);
        }
    }
}
