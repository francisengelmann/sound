// OpenAL headers
#include <al.h>
#include <alc.h>

// C/C++ includes
#include <list>
#include <iostream>
#include <vector>
#include <string>

// OpenCV includes
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

// FFTW includes
#include <fftw3.h>

/**
 * @brief get_audio_devices
 * @return std::vector containing name as std::string of input devices
 */
std::vector<std::string> get_audio_devices() {
  const ALCchar *devices = alcGetString(NULL, ALC_DEVICE_SPECIFIER);
  const ALCchar *device = devices, *next = devices + 1;
  size_t len = 0;
  std::vector<std::string> dev;

  while (device && *device != '\0' && next && *next != '\0') {
    dev.push_back(device);
    len = strlen(device);
    device += (len + 1);
    next += (len + 2);
  }

  return dev;
}

bool is_enumerate_supported() {
  ALboolean enumeration;
  enumeration = alcIsExtensionPresent(NULL, "ALC_ENUMERATION_EXT");
  return !(enumeration == AL_FALSE);
}

/*void plot_spectrum(double* buffer[2], int nsamples) {

  // Compute maximum amplitude
  bool scaling = true;
  int max = 0;
  if (scaling) {
    for (int i=0; i<nsamples; i++)
      if (max < buffer[i][0]) max = buffer[i][0];
    max*=2.0;
  } else {
    max = 65536/2;
  }

  // Plot points
  cv::Mat image(600,nsamples, CV_8UC3, cv::Scalar(0,0,0));
  std::vector<cv::Point> contour_real, contour_comp;
  contour_real.resize(nsamples+1);
  contour_comp.resize(nsamples+1);
  for (int i=0; i<nsamples; i++) {
    int it = (i+nsamples/2)%nsamples;
    int v_real = (int)((double)buffer[it][0]*600.0f/max)+300;
    int v_comp = (int)((double)buffer[it][0]*600.0f/max)+300;
    contour_real.push_back( cv::Point(i, v_comp) );
    contour_comp.push_back( cv::Point(i, v_real) );
  }

  const cv::Point *pts_real = (const cv::Point*) cv::Mat(contour_real).data;
  const cv::Point *pts_comp = (const cv::Point*) cv::Mat(contour_comp).data;
  int npts_real = cv::Mat(contour_real).rows;
  int npts_comp = cv::Mat(contour_comp).rows;

  cv::polylines(image, &pts_real, &npts_real, 1,
              false, 			// draw closed contour (i.e. joint end to start)
              cv::Scalar(0,255,0),// colour RGB ordering (here = green)
              1, 		        // line thickness
              CV_AA, 0);
  cv::polylines(image, &pts_comp, &npts_comp, 1,
              false, 			// draw closed contour (i.e. joint end to start)
              cv::Scalar(0,255,0),// colour RGB ordering (here = green)
              1, 		        // line thickness
              CV_AA, 0);

  cv::imshow("spectrum", image);
  cv::waitKey(1);
}*/

void plot_samples(short* buffer, int nsamples) {

  // Compute maximum amplitude
  bool scaling = false;
  int max = 0;
  if (scaling) {
    for (int i=0; i<nsamples; i++)
      if (max < buffer[i]) max = buffer[i];
    max*=2.0;
  } else {
    max = 65536/2;
  }

  // Plot points
  cv::Mat image(600,nsamples, CV_8UC3, cv::Scalar(0,0,0));
  std::vector<cv::Point> contour;
  contour.resize(nsamples+1);
  contour.clear();
  for (int i=1; i<nsamples; i++) {
    //int v = (int)((double)buffer[i]*600.0f/65536.0)+300;
    int v = (int)((double)buffer[i]*600.0f/max)+300;
    contour.push_back( cv::Point(i, v) );
  }

  const cv::Point *pts = (const cv::Point*) cv::Mat(contour).data;
  int npts = cv::Mat(contour).rows;

  cv::polylines(image, &pts, &npts, 1,
              false, 			// draw closed contour (i.e. joint end to start)
              cv::Scalar(0,255,0),// colour RGB ordering (here = green)
              2, 		        // line thickness
              CV_AA, 0);

  cv::imshow("spectrum", image);
  cv::waitKey(1);
}

int main(int argc, char ** argv) {

  // Check if we can list microphones and print the list
  std::vector<std::string> devs;
  if (is_enumerate_supported() ) {
    devs = get_audio_devices();
    for (auto s : devs) std::cout << s << std::endl;
  } else {
    std::cout << "Enumeration is not supported." << std::endl;
    return 0;
  }
  alGetError();

  // Open listening device
  ALCuint frequency = 44100/4;
  ALCenum format = AL_FORMAT_MONO16;
  ALCsizei buffer_size = 1024; // Anzahl der Sampleframes
  ALCchar* deviceName = NULL;
  ALCdevice *device = alcCaptureOpenDevice(deviceName, frequency, format, buffer_size);
  if (alGetError() != AL_NO_ERROR) return 0;

  // Set FFT
  int N = buffer_size;
  fftw_complex *in, *out;
  fftw_plan p;
  in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
  p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE_PATIENT);

  // Start capturing
  alcCaptureStart(device);

  int nsamples = buffer_size; // number of samples to take in one go, same as size of buffer
  int back = 3; // number of buffers to use for computing average
  double waterfall[back][nsamples/2]; // buffer for fourrier output, only have the size, first have of spectrum is symetric to second part
  double average[nsamples/2];

  int j = 0;
  while (true) {
    j++;
    j=j%back;

    // Check how many samples are available
    ALCenum param = ALC_CAPTURE_SAMPLES; // Gibt die Anzahl der verfügbaren Aufnahmesamples zurück
    ALCsizei size = (ALCsizei)sizeof(ALint); // Größe des bereitgestellten Puffers, hier ein int
    ALint available_samples; //hier wird der wert zureuckgeliefert

    alcGetIntegerv(device, param, size, &available_samples);
    if (alGetError() != AL_NO_ERROR) return 0;

    // Put the available samples into the buffer
    ALshort buffer[buffer_size];
    alcCaptureSamples(device, (ALCvoid *)buffer, available_samples);

    // Put buffer into format of fft
    for (int i=0; i<buffer_size; i++) {
      in[i][0] = (double)buffer[i]; // Real part
      in[i][1] = (double)0.0; // Imaginary part
    }

    // Perform the FFT
    fftw_execute(p); /* repeat as needed */ // exactly, so the above stuff should only be executed once!!

    // Plot the buffer
    plot_samples(buffer, buffer_size);

    // PLOT SPECTRUM
    {

      double v; // max freq value
      for (int i=0; i<nsamples/2; i++){
        v = std::sqrt(out[i][0]*out[i][0] + out[i][1]*out[i][1]);
        waterfall[j][i] = v;
      }

      bool scaling = false;
      int max = 0;
      int max_freq;
      for (int i=0; i<nsamples/2; i++){
        double avg = 0;
        for (int k=0; k<back; k++) {
          avg += waterfall[k][i];
        }
        average[i] = avg/back;
        if (max < average[i]) {
          max = average[i];
          max_freq = i;
        }
      }
      max = 65536*32;

      // Compute maximum amplitude

      //if (scaling) {
      /*  for (int i=0; i<nsamples/2; i++){
          v = std::sqrt(out[i][0]*out[i][0] + out[i][1]*out[i][1]);
          if (max < v) {
            max = v;
            max_freq = i;
          }
        }
        //max_freq=max;
      //} else {
      //}*/

      double cut_off_factor = 4;//

      // Plot points
      cv::Mat image(600,nsamples, CV_8UC3, cv::Scalar(20,30,40));

      for (int i=0; i<buffer_size; i+=40) cv::line(image, cv::Point(i,0), cv::Point(i,600), cv::Scalar(100, 100, 100), 1);
      for (int i=0; i<600; i+=40) cv::line(image, cv::Point(0,i), cv::Point(buffer_size,i), cv::Scalar(100, 100, 100), 1);

      std::vector<cv::Point> contour;
      contour.resize(nsamples+1);
      contour.clear();

      std::vector<int> maximums;
      maximums.clear();

      std::vector<cv::Point> contour_deriv;
      contour_deriv.resize(nsamples+1);
      contour_deriv.clear();
      contour_deriv.push_back(cv::Point(0,300));
      for (int i=0; i<nsamples/cut_off_factor; i++) {
        int it = i;
        double mag = average[i];
        double pos_y = (int)(mag*-598.0f/max)+598;
        contour.push_back( cv::Point(cut_off_factor*i, pos_y) );
        if (i>1) {

          double deriv_curr = average[i] - average[i-1];
          double deriv_prev = average[i-1] - average[i-2];

          if (std::signbit(deriv_curr)==1 && std::signbit(deriv_prev)==0 && deriv_prev - deriv_curr > 40000) {
            maximums.push_back(i-1);
            //std::cout << average[i] << " " << average[i-1] << "; \t";
            //std::cout << std::signbit(deriv_prev) << " " << deriv_prev << " - ";
            //std::cout << std::signbit(deriv_curr) << " " << deriv_curr << std::endl;
          }

          double pos_y = (int)(deriv_curr*-598.0f/max)+300;
          contour_deriv.push_back( cv::Point(cut_off_factor*i, pos_y ));
        }
      }

      // Draw max peak line with frequency
      /*if (v>17) {
        cv::line(image, cv::Point(cut_off_factor*max_freq,0), cv::Point(cut_off_factor*max_freq,600), cv::Scalar(255,255,255),2);
        cv::putText(image, std::to_string((max_freq*frequency)/buffer_size)+" Hz", cv::Point(50,50), 1, 2, cv::Scalar(255,255,255),1);
      }*/

      std::cout << "---------------" << std::endl;
      int i_min = 0;
      int fundamental=0;
      for (int i=0; i<maximums.size(); i++) {
        if (maximums.at(i) < maximums.at(i_min)) {
          i_min = i;
          fundamental = (maximums[i_min]*frequency)/buffer_size;
        }
      }
      if (maximums.size() > 0)
      std::cout << i_min << " " << maximums.at(i_min) << std::endl;

      cv::putText(image, std::to_string(fundamental)+" Hz", cv::Point(50,50), 1, 2, cv::Scalar(255,255,255),1);


      /*if (maximums.size() > 1) {
        double fundamental = 0.0;
        for (int i=1; i<maximums.size(); i++) {

          double diff = maximums.at(i)-maximums.at(i-1);
          std::cout << diff << std::endl;
          fundamental += diff;
        }
        fundamental /= maximums.size();
        fundamental = (fundamental*frequency)/buffer_size;
        cv::putText(image, std::to_string(fundamental)+" Hz", cv::Point(50,50), 1, 2, cv::Scalar(255,255,255),1);
      }*/

      for (auto i : maximums) {
        cv::line(image, cv::Point(cut_off_factor*i,0), cv::Point(cut_off_factor*i,600), cv::Scalar(255,200,200),2);
      }

      // Draw spectrum
      const cv::Point *pts = (const cv::Point*) cv::Mat(contour).data;
      int npts = cv::Mat(contour).rows;
      cv::polylines(image, &pts, &npts, 1,
                  false,             // Draw closed contour (i.e. joint end to start)
                  cv::Scalar(255,100,50), // Colour RGB ordering (here = green)
                  2,                 // Line thickness
                  0, 0);

      // Draw derivation
      const cv::Point *pts_deriv = (const cv::Point*) cv::Mat(contour_deriv).data;
      int npts_deriv = cv::Mat(contour_deriv).rows;
      cv::polylines(image, &pts_deriv, &npts_deriv, 1,
                  false,             // Draw closed contour (i.e. joint end to start)
                  cv::Scalar(50,100,250), // Colour RGB ordering (here = green)
                  2,                 // Line thickness
                  0, 0);

      cv::imshow("Spectrum", image);
      cv::waitKey(0);
    }
  }

  fftw_destroy_plan(p);
  fftw_free(in); fftw_free(out);

  alcCaptureStop(device);
  alcCaptureCloseDevice(device);
}

/*int main2(int argC,char* argV[])
{
  std::list<ALuint> bufferQueue; // A quick and dirty queue of buffer objects

  ALenum errorCode=0;
  ALuint helloBuffer[16], helloSource[1];
  ALCdevice* audioDevice = alcOpenDevice(NULL); // Request default audio device
  errorCode = alcGetError(audioDevice);
  ALCcontext* audioContext = alcCreateContext(audioDevice,NULL); // Create the audio context
  alcMakeContextCurrent(audioContext);
  errorCode = alcGetError(audioDevice);
  
  // Request the default capture device with a half-second buffer
  ALCdevice* inputDevice = alcCaptureOpenDevice(NULL,FREQ,AL_FORMAT_MONO16,FREQ/2);
  errorCode = alcGetError(inputDevice);
  alcCaptureStart(inputDevice); // Begin capturing
  errorCode = alcGetError(inputDevice);

  alGenBuffers(16,&helloBuffer[0]); // Create some buffer-objects
  errorCode = alGetError();

  // Queue our buffers onto an STL list
  for (int ii=0;ii<16;++ii) {
    std::cout << helloBuffer[ii] << std::endl;
    bufferQueue.push_back(helloBuffer[ii]);
  }


  alGenSources (1, &helloSource[0]); // Create a sound source
  errorCode = alGetError();
}*/
