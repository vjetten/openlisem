/*************************************************************************
**  openLISEM: a spatial surface water balance and soil erosion model
**  Copyright (C) 2010,2011, 2020  Victor Jetten
**  contact: v.g.jetten AD utwente DOT nl
**
**  This program is free software: you can redistribute it and/or modify
**  it under the terms of the GNU General Public License GPLv3 as published by
**  the Free Software Foundation, either version 3 of the License, or
**  (at your option) any later version.
**
**  This program is distributed in the hope that it will be useful,
**  but WITHOUT ANY WARRANTY; without even the implied warranty of
**  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
**  GNU General Public License v3 for more details.
**
**  You should have received a copy of the GNU General Public License GPLv3
**  along with this program.  If not, see <http://www.gnu.org/licenses/>.
**
**  Authors: Victor Jetten, Bastian van de Bout
**  Developed in: MingW/Qt/
**  website, information and code: http://lisem.sourceforge.net
**
*************************************************************************/


#pragma once
#include <array>
#include <vector>


template<
    typename T,
    size_t nr_dimensions>
class Array
{

public:

                   Array               ();

                   Array               (std::initializer_list<T> elements);

                   Array               (size_t size);

                   Array               (size_t size1,
                                        size_t size2);

                   Array               (Array const& other)=default;

                   Array               (Array&& other)=default;

    virtual        ~Array              ()=default;

    Array&         operator=           (Array const& other)=default;

    Array&         operator=           (Array&& other)=default;

    size_t         size                () const;

    T*             data                ();

    T&             operator[]          (size_t index);

    T const&       operator[]          (size_t index) const;

    std::array<size_t, nr_dimensions> const&
                   shape               () const;

private:

    std::vector<T> _cells;

    std::array<size_t, nr_dimensions> _shape;

};


template<
    typename T,
    size_t nr_dimensions>
inline Array<T, nr_dimensions>::Array()

    : _cells(),
      _shape()

{
}


template<
    typename T,
    size_t nr_dimensions>
inline Array<T, nr_dimensions>::Array(
    std::initializer_list<T> elements)

    : _cells(elements),
      _shape{{elements.size()}}

{
    static_assert(nr_dimensions == 1, "");
}


template<
    typename T,
    size_t nr_dimensions>
inline Array<T, nr_dimensions>::Array(
    size_t size)

    : _cells(size),
      _shape{{size}}

{
    static_assert(nr_dimensions == 1, "");
}


template<
    typename T,
    size_t nr_dimensions>
inline Array<T, nr_dimensions>::Array(
    size_t size1,
    size_t size2)

    : _cells(size1 * size2),
      _shape{{size1, size2}}

{
    static_assert(nr_dimensions == 2, "");
}


template<
    typename T,
    size_t nr_dimensions>
inline size_t Array<T, nr_dimensions>::size() const
{
    return _cells.size();
}


template<
    typename T,
    size_t nr_dimensions>
inline T* Array<T, nr_dimensions>::data()
{
    return _cells.data();
}


template<
    typename T,
    size_t nr_dimensions>
inline T& Array<T, nr_dimensions>::operator[](
        size_t index)
{
    return _cells[index];
}


template<
    typename T,
    size_t nr_dimensions>
inline T const& Array<T, nr_dimensions>::operator[](
        size_t index) const
{
    return _cells[index];
}


template<
    typename T,
    size_t nr_dimensions>
inline std::array<size_t, nr_dimensions> const&
        Array<T, nr_dimensions>::shape() const
{
    return _shape;
}
