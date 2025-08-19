# Golf Launch Monitor - Radar DSP Algorithm Development

## Milestone 1: Data Foundation & Analysis

This MATLAB project implements the first milestone of the golf launch monitor radar DSP algorithm development for improved club/ball separation accuracy.

### Project Overview

- **Objective**: Develop advanced radar signal processing algorithms for golf launch monitor
- **Target**: ±1 mph accuracy vs TrackMan for both club and ball speeds
- **Focus**: Improved wedge shot accuracy and false positive elimination
- **Hardware**: InnoSenT SMR333 24 GHz CW Doppler radar module

### Milestone 1 Deliverables

**Acceptance Criteria:**
- ✅ Dataset fully parsed (≥90% of test shots)
- ✅ FFT/spectrogram alignment with TrackMan timing (≥90%)
- ✅ Baseline STFT processing pipeline established

### File Structure

```
Golf_Launch_Monitor_DSP/
├── main_analysis.m                    # Main analysis script
├── parse_golf_dataset.m               # Dataset parsing functions
├── load_iq_data.m                     # I/Q data loader
├── baseline_stft_processing.m         # STFT processing pipeline
├── select_representative_shots.m      # Shot selection utilities
├── validate_against_trackman.m        # TrackMan validation framework
├── analyze_signal_quality.m           # Signal quality assessment
├── create_analysis_visualizations.m   # Visualization generation
├── generate_milestone1_report.m       # Report generation
├── README.md                          # This file
├── milestone1_results.mat             # Analysis results (generated)
├── milestone1_report.txt              # Completion report (generated)
└── milestone1_figures/                # Visualizations (generated)
    ├── spectrograms.png
    ├── signal_quality.png
    ├── trackman_validation.png
    ├── velocity_detection.png
    └── data_quality.png
```

### Usage Instructions

1. **Set Data Path**: Ensure your dataset is in the correct location:
   ```matlab
   config.data_path = 'D:\MATLAB OR PYTHON  PROJECT\28-01-2024 DATA';
   ```

2. **Run Analysis**: Execute the main script:
   ```matlab
   run('main_analysis.m');
   ```

3. **Review Results**: Check generated files:
   - `milestone1_results.mat` - Complete analysis results
   - `milestone1_report.txt` - Acceptance criteria evaluation
   - `milestone1_figures/` - Analysis visualizations

### Key Features

**Data Processing:**
- JSON metadata parsing
- 4-channel I/Q binary data loading
- Comprehensive error handling and validation

**Signal Processing:**
- STFT with Hamming windowing
- Adaptive peak detection
- Multi-target tracking framework
- Velocity estimation using Doppler conversion

**Validation Framework:**
- TrackMan reference comparison
- Timing alignment assessment
- Preliminary accuracy metrics
- Category-specific analysis (wedge/iron/driver)

**Quality Assessment:**
- SNR analysis
- Noise floor estimation
- Saturation detection
- DC offset monitoring

### Configuration Parameters

```matlab
config.sampling_freq = 22700;           % Hz
config.speed_coef = 0.1388888888888889; % Frequency to velocity conversion
config.stft_window = 256;               % STFT window size
config.stft_overlap = 128;              % STFT overlap
config.fft_size = 1024;                 % FFT points
config.channels = 4;                    % I/Q channels
```

### Expected Results

**Milestone 1 Targets:**
- Dataset parsing success: ≥90%
- TrackMan timing alignment: ≥90%
- Signal quality score: ≥0.7
- Peak detection rate: ≥80%

### Next Steps (Milestone 2)

1. **Adaptive CFAR Detection**: Implement robust target detection
2. **Multi-Target Tracking**: Develop sophisticated tracking algorithms
3. **Club/Ball Separation**: Core differentiation logic
4. **Wedge Optimization**: Enhanced processing for low-speed shots

### Technical Notes

**Radar Specifications:**
- Module: InnoSenT SMR333
- Frequency: 24 GHz CW Doppler
- Sampling: 22.7 kHz
- Data format: 16-bit signed integers

**MCU Considerations:**
- Current: Renesas Rx65 (2MB Flash, 640KB SRAM)
- Future: STM32U (512KB Flash, 274KB SRAM)
- Algorithm designed for eventual embedded deployment

### Contact

For questions or issues regarding this milestone implementation, please refer to the project documentation or contact the development team.

---

**Project**: Golf Launch Monitor Club/Ball Separation & Accuracy Improvement  
**Developer**: Tanveer Hussain  
**Date**: January 2025  
**Version**: 1.0
