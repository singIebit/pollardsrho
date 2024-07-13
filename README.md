# Pollard's Rho and Brent's Algorithm for SECP256K1 Curve

![C](https://img.shields.io/badge/language-C-blue)
![GMP](https://img.shields.io/badge/dependency-GMP-green)
![pthreads](https://img.shields.io/badge/dependency-pthreads-yellow)

## Description

Pollard's rho algorithm is renowned for integer factorization, leveraging minimal space with an expected runtime proportional to the square root of the smallest prime factor. Brent's cycle-detection algorithm, originating from numerical analysis, enhances Pollard's rho by efficiently identifying cycles in sequences of iterates. This hybrid method combines the bisection, secant, and inverse quadratic interpolation techniques.

This repository features implementations that optimize Pollard's rho with Brent's algorithm to tackle discrete logarithm challenges on the secp256k1 elliptic curve. The goal is to determine the relationship between points 𝐺 (base point),  𝑃 (given point),  and the integer 𝑥 such that 𝑃 = 𝑥𝐺.

## Features

- **Collision Detection**
  - Finds collisions efficiently using Pollard's Rho and Brent's Algorithm.
  - Optimized for the secp256k1 elliptic curve.

- **Parallel Processing**
  - Utilizes multithreading with pthreads to speed up calculations.

- **Progress Monitoring**
  - Displays real-time progress updates during the search for collisions.

## Installation

1. Clone this repository:
    ```bash
    git clone https://github.com/singIebit/pollardsrho.git
    cd pollardsrho
    ```

2. Install the libgmp-dev:
    ```bash
    sudo apt-get update
    sudo apt-get install libgmp-dev
    ```

3. Compile the project:
    ```bash
    cmake .
    make
    ```

4. Run the program:
    ```bash
    ./pollardsrho <public key> <key range>
    ```

    Replace `<public key>` with the point \( G \) on the secp256k1 curve multiplied by your private key value, and `<key range>` with the size of the search interval for \( k \).

    Example usage:
    ```bash
    ./pollardsrho 03633cbe3ec02b9401c5effa144c5b4d22f87940259634858fc7e59b1c09937852 130
    ```

## Commands

- `./pollardsrho <public key> <key range>`: Starts the search for the private key corresponding to the given public key within the specified range.

## Technologies Used

- C Programming Language
- GMP (GNU Multiple Precision Arithmetic Library)
- pthreads (POSIX Threads Library)

## Contributing

Contributions are welcome! Feel free to open issues or submit pull requests.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

[![Add Star](https://img.shields.io/github/stars/singIebit/pollardsrho?style=social)](https://github.com/singIebit/pollardsrho)

---

<p align="center">
  <a href="https://github.com/singIebit">
    <img src="https://readme-typing-svg.demolab.com?font=Georgia&size=18&duration=2000&pause=100&multiline=true&width=500&height=80&lines=No+Clip+Studio;Programmer+%7C+Student+%7C+Cyber+Security;+%7C+Android+%7C+Apps" alt="Typing SVG" />
  </a>
</p>

<a href="https://github.com/singIebit">
    <img src="https://github-stats-alpha.vercel.app/api?username=singIebit&cc=22272e&tc=37BCF6&ic=fff&bc=0000">
</a>

- 🔭 I’m currently working on [Data Leak Search](https://play.google.com/store/apps/details?id=com.NoClipStudio.DataBaseSearch)

- 🚀 I’m looking to collaborate on: [Cyber-Security](https://play.google.com/store/apps/details?id=com.hashsuite.droid)

- 📝 I regularly read: [https://nvd.nist.gov/vuln](https://nvd.nist.gov/vuln)

- 📄 Know about my experiences: [https://www.linkedin.com/in/lucas-jhonatan-215594208](https://www.linkedin.com/in/lucas-jhonatan-215594208)

<br>
My Github Stats

![](http://github-profile-summary-cards.vercel.app/api/cards/profile-details?username=singIebit&theme=dracula) 
![](http://github-profile-summary-cards.vercel.app/api/cards/repos-per-language?username=singIebit&theme=dracula) 
![](http://github-profile-summary-cards.vercel.app/api/cards/most-commit-language?username=singIebit&theme=dracula)

<h3 align="left">Connect with me:</h3>
<p align="left">
<a href="https://linkedin.com/in/lucas-jhonatan-215594208" target="blank"><img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/linked-in-alt.svg" alt="lucas-jhonatan-215594208" height="30" width="40" /></a>
<a href="https://www.youtube.com/@noclipstudiobr" target="blank"><img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/youtube.svg" alt="@noclipstudiobr" height="30" width="40" /></a>
<a href="https://discord.gg/https://discord.gg/wXqcJDHht8" target="blank"><img align="center" src="https://raw.githubusercontent.com/rahuldkjain/github-profile-readme-generator/master/src/images/icons/Social/discord.svg" alt="https://discord.gg/wXqcJDHht8" height="30" width="40" /></a>
</p>

<h3 align="left">Languages and Tools:</h3>
<p align="left"> <a href="https://developer.android.com" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/android/android-original-wordmark.svg" alt="android" width="40" height="40"/> </a> <a href="https://www.w3schools.com/cs/" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/csharp/csharp-original.svg" alt="csharp" width="40" height="40"/> </a> <a href="https://git-scm.com/" target="_blank" rel="noreferrer"> <img src="https://www.vectorlogo.zone/logos/git-scm/git-scm-icon.svg" alt="git" width="40" height="40"/> </a> <a href="https://www.linux.org/" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/linux/linux-original.svg" alt="linux" width="40" height="40"/> </a> <a href="https://www.python.org" target="_blank" rel="noreferrer"> <img src="https://raw.githubusercontent.com/devicons/devicon/master/icons/python/python-original.svg" alt="python" width="40" height="40"/> </a> <a href="https://unity.com/" target="_blank" rel="noreferrer"> <img src="https://www.vectorlogo.zone/logos/unity3d/unity3d-icon.svg" alt="unity" width="40" height="40"/> </a> </p>
