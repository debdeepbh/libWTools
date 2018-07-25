#pragma once

#include <fftw3.h>
#include <complex>
#include <valarray>

using std::complex;

// namespace FWTransform {

    ///
    /// \brief A time-domain signal
    ///
    using Signal = std::valarray<double>;

    ///
    /// \brief A frequency-domain representation of a signal.
    ///
    using Spectrum = std::valarray<complex<double>>;

    ///
    /// \brief
    ///
    class WTools{

    private:

        const int length; ///< The length of the time-data.
        // fftw_complex* in; ///< The memory to store the input data for FFTW execution.
        // fftw_complex* out; ///< The memory to store the output data for FFTW execution.
        // fftw_plan forward_plan; ///< The plan to execute forward FFTs
        // fftw_plan backward_plan; ///< The plan to execute forward FFTs

    public:

        ///
        /// \brief Construct a FWTransform object and associated FFTW plans
        ///
        /// All signals given to this object must be of length `length`, otherwise
        /// an exception will be thrown.
        ///
        WTools(const int _length) : length(_length) {};

            // we create one plan to force FFTW to cache plans of this length
            // fftw_
            //     fftw_plan caching_plan{fftw_plan_dft_1d(length, )}


                                         // in(reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*length))),
                                         // out(reinterpret_cast<fftw_complex*>(fftw_malloc(sizeof(fftw_complex)*length))),
                                         // forward_plan(fftw_plan_dft_1d(length, in, out, FFTW_FORWARD, FFTW_MEASURE)),
                                         // backward_plan(fftw_plan_dft_1d(length, in, out, FFTW_BACKWARD, FFTW_MEASURE)) {};
        ///
        /// \brief FWTransform destructor. Free plan objects and assigned memory.
        ///
        ~WTools() {

            // free all FFTW objects
            // fftw_destroy_plan(this->forward_plan);
            // fftw_destroy_plan(this->backward_plan);
            // fftw_free(in);
            // fftw_free(out);

        }

        ///
        /// \brief Convolve time-domain signals `A` and `B`.
        ///
        auto convolve(const Signal A, const Signal B) -> Signal;

        ///
        /// \brief A wrapper to perform convolution on raw double points
        ///
        /// This function does not check that the memory pointed to by `A` or `B`
        /// is of the right length. Passing `A` or `B` with the wrong length is
        /// undefined behaviour (UB).
        ///
        /// USER MUST FREE RETURNED MEMORY WHEN DONE.
        ///
        auto convolve(const double* A, const double* B) -> double*;

        ///
        /// \brief Convolve time-domain signals `A` and `B`.
        ///
        auto deconvolve(const Signal input, const Signal response) -> Signal;

        ///
        /// \brief A wrapper to perform convolution on raw double points
        ///
        auto deconvolve(const double* input, const double* response) -> double*;

	void fft(int size, fftw_complex* input, fftw_complex* output);


    }; // END: class FWTransform

//} // END: namespace FWTransform
