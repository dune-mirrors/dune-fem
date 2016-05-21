/*
    pybind11/reference_wrapper.h: Analogue to std::reference_wrapper for Python objects

    Copyright (c) 2016 Martin Nolte

    All rights reserved. Use of this source code is governed by a
    BSD-style license that can be found in the LICENSE file.
*/

#pragma once

#include <functional>
#include <type_traits>
#include <utility>

NAMESPACE_BEGIN(pybind11)

template <typename T> class reference_wrapper {
    template <typename U>
    friend class reference_wrapper;
public:
    reference_wrapper(object obj) : obj(std::move(obj)), ref(this->obj.cast<T&>()) {}
    reference_wrapper(T &t) : obj(detail::get_object_handle(&t), true), ref(&t) {}
    reference_wrapper(T &&) = delete;
    reference_wrapper(const reference_wrapper &) = default;
    reference_wrapper(reference_wrapper &&) = default;
    template <typename U> reference_wrapper(const reference_wrapper<U> &other) : obj(other.obj), ref(other.ref) {}
    template <typename U> reference_wrapper(reference_wrapper<U> &&other) : obj(std::move(other.obj)), ref(other.ref) {}
    reference_wrapper &operator=(object obj) {
        this->obj = std::move(obj);
        this->ref = &this->obj.cast<T&>();
        return *this;
    }
    reference_wrapper &operator=(T &t) {
        this->obj = object(detail::get_object_handle(&t), true);
        this->ref = &t;
        return *this;
    }
    reference_wrapper &operator=(const reference_wrapper &) = default;
    reference_wrapper &operator=(reference_wrapper &&) = default;
    template <typename U> reference_wrapper &operator=(const reference_wrapper<U> &other) {
        this->obj = other.obj;
        this->ref = other.ref;
        return *this;
    }
    template <typename U> reference_wrapper &operator=(reference_wrapper<U> &&other) {
        this->obj = std::move(other.obj);
        this->ref = other.ref;
        return *this;
    }
    operator object() const { return obj; }
    operator T&() const { return *ref; }
    operator std::reference_wrapper<T>() const { return std::reference_wrapper<T>(*ref); }
    T &get() const { return *ref; }
    template <typename... Args> typename std::result_of<T&(Args&&...)>::type operator()(Args&&... args) {
        return get()(std::forward<Args>(args)...);
    }
private:
    object obj;
    T *ref;
};

template <typename T> reference_wrapper<T> ref(T &t) { return reference_wrapper<T>(t); }
template <typename T> reference_wrapper<T> ref(reference_wrapper<T> t) { return std::move(t); }
template <typename T> reference_wrapper<const T> cref(T &t) { return reference_wrapper<const T>(t); }
template <typename T> reference_wrapper<const T> cref(reference_wrapper<T> t) { return std::move(t); }

NAMESPACE_END(pybind11)
