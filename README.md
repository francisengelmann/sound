## Sound

This little application samples and visualizes the input from the sound card.
Additionally, it shows the spectrum of the incoming signal and displays the frequency of the maximum span.

### Screenshots

![alt tag](data/screen1.png)
![alt tag](data/screen2.png)
Spectrum resulting from playing a guitar string, you can clearly see the the harmonics of the base tone!

### Dependencies
You need the following libraries:
- OpenCV to visualize the sound input and spectrum
- OpenAL to access the sound card
- FFTW to apply fourrier transformations

### Installation
Set up the correct pathes for the dependencies, then type this into your console:
```
mkdir build; cd build;
cmake ..; build -j;
```

#### References
[1] Daniel Shawcross Wilkerson. Harmony Explained: Progress Towards a Scientific Theory of Music
https://arxiv.org/html/1202.4212v1
