#pragma once
#include "array.h"


template<
    typename T>
class Raster:
    private Array<T, 2>
{

public:

                   Raster              ();

                   Raster              (size_t nr_rows,
                                        size_t nr_cols);

                   Raster              (Raster const& other)=default;

                   Raster              (Raster&& other)=default;

    virtual        ~Raster             ()=default;

    Raster&        operator=           (Raster const& other)=default;

    Raster&        operator=           (Raster&& other)=default;

    size_t         nr_rows             () const;

    size_t         nr_cols             () const;

    size_t         nr_cells            () const;

    size_t         index               (size_t row,
                                        size_t col) const;

    T*             operator[]          (size_t row);

    T const*       operator[]          (size_t row) const;

    T&             cell                (size_t index);

    T const&       cell                (size_t index) const;

    T&             cell                (size_t row,
                                        size_t col);

    T const&       cell                (size_t row,
                                        size_t col) const;

private:

};


template<
    typename T>
inline Raster<T>::Raster()

    : Array<T, 2>()

{
}


template<
    typename T>
inline Raster<T>::Raster(
    size_t nr_rows,
    size_t nr_cols)

    : Array<T, 2>(nr_rows, nr_cols)

{
}


template<
    typename T>
inline size_t Raster<T>::nr_rows() const
{
    return this->shape()[0];
}


template<
    typename T>
inline size_t Raster<T>::nr_cols() const
{
    return this->shape()[1];
}


template<
    typename T>
inline size_t Raster<T>::nr_cells() const
{
    return this->size();
}


template<
    typename T>
size_t Raster<T>::index(
    size_t row,
    size_t col) const
{
    return row * nr_cols() + col;
}


template<
    typename T>
inline T* Raster<T>::operator[](
    size_t row)
{
    return &(Array<T, 2>::operator[](index(row, 0)));
}


template<
    typename T>
inline T const* Raster<T>::operator[](
    size_t row) const
{
    return &(Array<T, 2>::operator[](index(row, 0)));
}


template<
    typename T>
inline T& Raster<T>::cell(
    size_t index)
{
    return Array<T, 2>::operator[](index);
}


template<
    typename T>
inline T const& Raster<T>::cell(
    size_t index) const
{
    return Array<T, 2>::operator[](index);
}


template<
    typename T>
inline T& Raster<T>::cell(
    size_t row,
    size_t col)
{
    return Array<T, 2>::operator[](index(row, col));
}


template<
    typename T>
inline T const& Raster<T>::cell(
    size_t row,
    size_t col) const
{
    return Array<T, 2>::operator[](index(row, col));
}
