# Reference-Anchored Envelope Modelling for Drought-Driven Tree Mortality Risk

**A reproducible Google Earth Engine framework that predicts tree-mortality hotspots in Mediterranean and semi-arid forests, at landscape to national scale.**

Dr. Filippos Eliades · Remote Sensing and GeoEnvironment Lab, Cyprus University of Technology · Eratosthenes Centre of Excellence · 2026

---

## What it does

This framework flags where drought-driven tree mortality is likely across a forest domain. Instead of assessing trees individually, it transfers a mortality "signature" from a small set of confirmed mortality footprints to the wider landscape, producing an actionable high-risk screening layer.

It was calibrated on three confirmed mortality sites in Cyprus (Akamas, Stavrovouni, Machairas), spanning contrasting species, terrain and substrate, and anchored to the 2008 and 2016 drought-mortality events.

## Prospective validation

In 2025, two independent mortality sites appeared in Cyprus. Both fell inside zones the framework had already mapped as high-risk using only the 2008/2016 signatures, with 83.5% and 82.7% of their forested area exceeding the high-intensity threshold. With 18.88% of Cyprus's forest flagged as high-risk, the probability of this co-occurrence being random was p ≈ 0.036 (binomial test). This is a genuine out-of-sample confirmation, not a retrospective fit.

## Example output

<img width="100%" height="5100" alt="Layout4" src="https://github.com/user-attachments/assets/5c38b766-0d2e-425e-988c-1b1b049ce135" />

## How it works

The method combines four components:

- **Physiographic similarity** — quantile envelopes (0.05–0.95) for slope, multi-scale topographic position (TPI 300 m / 1000 m) and curvature, combined under a majority rule (≥60% of criteria)
- **Exposure constraints** — circular aspect windows, warm-season solar loading (April–September), and lithology filtering
- **Drought gate** — SPEI-12 ≤ −1.4
- **Spectral collapse** — Landsat event-year decline (NBR / NDVI) against a fixed 3-year pre-event baseline, with a year+1 persistence rule to suppress transient stress

Outputs are STRICT/RELAX resemblance masks and a kernel-density HEAT intensity surface (0–1) for polygon-level risk scoring.

## Data and platform

- **Platform:** Google Earth Engine (JavaScript API), run at 30 m
- **Inputs (all public via GEE):** Landsat Collection 2 Level-2 (USGS), Copernicus DEM GLO-30 (ESA), CSIC SPEI v2.10, national forest mask, geological substrate

## Published research

- Eliades et al. (2026). *Reference-anchored envelope modelling for predicting tree mortality hotspots driven by drought.* International Journal of Applied Earth Observation and Geoinformation 154, 105609. https://doi.org/10.1016/j.jag.2026.105609
- Eliades et al. (2026). *Forests in a semi-arid climate die with a memory.* Journal of Forestry Research 37, 79. https://doi.org/10.1007/s11676-026-02016-z
- Eliades et al. (2024). *Understanding Tree Mortality Patterns.* Forests 15, 1357. https://doi.org/10.3390/f15081357

## Citation

If you use this framework, please cite the IJAEOG (2026) paper above and the archived code:

> Eliades, F. (2026). Zenodo. https://doi.org/10.5281/zenodo.20154761

## License

Released under **CC BY-NC 4.0** — free to use and adapt for non-commercial purposes with attribution. For commercial use, please get in touch.

## Contact

Dr. Filippos Eliades
filippos.eliades@gmail.com · [LinkedIn](https://www.linkedin.com/in/filippos-%CE%B5liades-381627134/))

*Available for consultancy on forest mortality risk mapping in Mediterranean and semi-arid ecosystems.*
