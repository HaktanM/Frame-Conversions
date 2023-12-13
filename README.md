# Frame-Conversions
Earth Centered Earth Fixed (ECEF) Frame - Geodetic Frame Conversion Utils 

This repo aims to create user-friendly cpp scripts for frame conversions. 
Convert latitude-longitute-altitude (geodetic coordinates) to ECEF coordinates or vice versa. 

Make sure that you have installed CMake. You can easily download the latest version from
[web page](https://cmake.org/download/).


If you are not familiar with cmake, you can follow the first two steps [here](https://cmake.org/cmake/help/latest/guide/tutorial/index.html#introduction).


Even though you are not familiar with cmake, you can run the code via following steps
```sh
cd your_repo      
cmake .            # Start a cmake environment
cmake --build .    # Build the code
./convert_frames   # Run the executable
