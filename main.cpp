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

void plot_samples(ALshort* buffer, int nsamples=22050) {

  int max = 0;
  for (int i=0; i<nsamples; i++) if (max < buffer[i]) max = buffer[i];

  cv::Mat image(600,nsamples, CV_8UC3, cv::Scalar(0,0,0));
  std::vector<cv::Point> contour;
  contour.resize(nsamples+1);
  for (int i=0; i<nsamples; i++) {
    //int v = (int)((double)buffer[i]*600.0f/65536.0)+300;
    int v = (int)((double)buffer[i]*600.0f/max/2)+300;
    contour.push_back( cv::Point(i, v) );
  }

  const cv::Point *pts = (const cv::Point*) cv::Mat(contour).data;
  int npts = cv::Mat(contour).rows;

  cv::polylines(image, &pts, &npts, 1,
              false, 			// draw closed contour (i.e. joint end to start)
              cv::Scalar(0,255,0),// colour RGB ordering (here = green)
              1, 		        // line thickness
              CV_AA, 0);

  cv::imshow("plot", image);
  cv::waitKey(10);
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
  ALCuint frequency = 44100;
  ALCenum format = AL_FORMAT_MONO16;
  ALCsizei buffer_size = 1024; // Anzahl der Sampleframes
  ALCchar* deviceName = NULL;
  ALCdevice *device = alcCaptureOpenDevice(deviceName, frequency, format, buffer_size);
  if (alGetError() != AL_NO_ERROR) return 0;

  // Start capturing
  alcCaptureStart(device);

  while (true) {

    // Check how many samples are available
    ALCenum param = ALC_CAPTURE_SAMPLES; // Gibt die Anzahl der verfügbaren Aufnahmesamples zurück
    ALCsizei size = (ALCsizei)sizeof(ALint); // Größe des bereitgestellten Puffers, hier ein int
    ALint available_samples; //hier wird der wert zureuckgeliefert
    alcGetIntegerv(device, param, size, &available_samples);
    if (alGetError() != AL_NO_ERROR) return 0;

    // Put the available samples into the buffer
    ALshort buffer[1024];
    alcCaptureSamples(device, (ALCvoid *)buffer, available_samples);

    // Plot the buffer
    plot_samples(&buffer[0], 1024);
  }

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
