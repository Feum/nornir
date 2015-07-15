/*
 * utils.hpp
 *
 * Created on: 09/07/2015
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of AdaptiveFastFlow.
 *
 *  AdaptiveFastFlow is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  AdaptiveFastFlow is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with AdaptiveFastFlow.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

/*!
 * \file utils.hpp
 * \brief Implementation of vairous utilities.
 **/

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include <vector>

namespace adpff{

template<typename T> class Window{
private:
    std::vector<T> _windowImpl;
    size_t _span;
    size_t _nextIndex;
    size_t _storedValues;
public:
    Window(size_t span):_span(span), _nextIndex(0), _storedValues(0){
        _windowImpl.resize(_span);
    }

    void add(const T& value){
        _windowImpl.at(_nextIndex) = value;
        _nextIndex = (_nextIndex + 1) % _span;
        if(_storedValues < _span){
            ++_storedValues;
        }
    }

    std::vector<T> toVector() const{
        return _windowImpl;
    }

    void reset(){
        _windowImpl.clear();
        _windowImpl.resize(_span);
        _nextIndex = 0;
        _storedValues = 0;
    }

    size_t size() const{
        return _storedValues;
    }

    T& operator[](size_t idx){
        return _windowImpl[idx];
    }

    const T& operator[](size_t idx) const{
        return _windowImpl[idx];
    }

    T& at(size_t idx){
        return _windowImpl.at(idx);
    }

    const T& at(size_t idx) const{
        return _windowImpl.at(idx);
    }

    T sum() const{
        T r;

        for(size_t i = 0; i < _storedValues; i++){
            r += _windowImpl[i];
        }

        return r;
    }

    T average() const{
        T r;
        if(!_storedValues){
            return r;
        }

        r = sum();
        return r / _storedValues;
    }

};

}

#endif /* UTILS_HPP_ */
