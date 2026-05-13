// 3 AOI
var area1FC = ee.FeatureCollection('users/filipposeliades31/Akamas');
var area2FC = ee.FeatureCollection('users/filipposeliades31/Stavrovouni');
var area3FC = ee.FeatureCollection('users/filipposeliades31/Machairas');

// Unified geometry per AOI
var gArea1 = area1FC.geometry();
var gArea2 = area2FC.geometry();
var gArea3 = area3FC.geometry();

// Time range for SPEI
var startYear = 1991;
var endYear   = 2022;

// =====================
// Forest mask
// =====================
var WE_forestsFC    = ee.FeatureCollection('users/filipposeliades31/CyprusForestsAll');
var WE_forestsGeom  = WE_forestsFC.geometry();
var WE_forestMask30 = ee.Image().byte().paint({featureCollection: WE_forestsFC, color: 1}).rename('forest');


/***** Landsat helpers (για NDVI/NBR & drop/persist) *****/

// Scale factors για Landsat C2 L2 (optical)
function applyScaleFactors(image) {
  var optical = image.select('SR_B.*').multiply(0.0000275).add(-0.2);
  return image.addBands(optical, null, true);
}

// Rename to common names
function renameOli(img) { // Landsat 8/9
  return img.select(
    ['SR_B2','SR_B3','SR_B4','SR_B5','SR_B6','SR_B7','QA_PIXEL'],
    ['blue','green','red','nir','swir1','swir2','QA_PIXEL']
  );
}
function renameEtm(img) { // Landsat 5/7
  return img.select(
    ['SR_B1','SR_B2','SR_B3','SR_B4','SR_B5','SR_B7','QA_PIXEL'],
    ['blue','green','red','nir','swir1','swir2','QA_PIXEL']
  );
}

// Cloud & shadow mask (QA_PIXEL: bit 3 cloud, bit 4 shadow)
function fmask(img) {
  var qa = img.select('QA_PIXEL');
  var noCloud  = qa.bitwiseAnd(1 << 3).eq(0);
  var noShadow = qa.bitwiseAnd(1 << 4).eq(0);
  return img.updateMask(noCloud.and(noShadow));
}

// Indicators
function addNDVI(img) {
  var ndvi = img.normalizedDifference(['nir','red']).rename('NDVI');
  return img.addBands(ndvi);
}
function addNBR(img) {
  var nbr = img.normalizedDifference(['nir','swir2']).rename('NBR');
  return img.addBands(nbr);
}

// Pipelines: production of NDVI & NBR
function prepOli(img) {
  return addNBR(addNDVI(fmask(renameOli(applyScaleFactors(img)))))
    .select(['NDVI','NBR'])
    .copyProperties(img, ['system:time_start']);
}
function prepEtm(img) {
  return addNBR(addNDVI(fmask(renameEtm(applyScaleFactors(img)))))
    .select(['NDVI','NBR'])
    .copyProperties(img, ['system:time_start']);
}

// Landsat merged helper
function WE_mergedLandsatYearly(start, end, roi){
  var filter = ee.Filter.and(ee.Filter.bounds(roi), ee.Filter.date(start, end));
  var L5 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filter(filter).map(prepEtm);
  var L7 = ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filter(filter).map(prepEtm);
  var L8 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filter(filter).map(prepOli);
  var L9 = ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filter(filter).map(prepOli);
  return L5.merge(L7).merge(L8).merge(L9);
}
function WE_yearMeanImage(band, year, roi){
  var s = ee.Date.fromYMD(year,1,1);
  var e = s.advance(1,'year');
  return WE_mergedLandsatYearly(s, e, roi).select(band).mean();
}

/***** NDVI / NBR ANNUAL MEAN CHARTS *****/

// Area1 → NBR
// Area2 → NBR
// Area3 → NDVI

// Helper: creates an annual-mean ImageCollection for a band
function yearlyMeanIC(band, roi, y0, y1){
  var years = ee.List.sequence(y0, y1);

  var images = years.map(function(y){
    y = ee.Number(y);
    var im = WE_yearMeanImage(band, y, roi)
      .rename(band)
      .set('year', y)
      .set('system:time_start', ee.Date.fromYMD(y, 7, 1).millis());
    return im;
  });

  return ee.ImageCollection(images);
}

// Generic helper για chart
function chartIndexSeries(ic, bandName, roi, title, yLabel, scale){
  scale = scale || 30;

  var chart = ui.Chart.image.series({
      imageCollection: ic.select(bandName),
      region: roi,
      reducer: ee.Reducer.mean(),
      scale: scale
    })
    .setChartType('LineChart')
    .setOptions({
      title: title,
      hAxis: { title: 'Year' },
      vAxis: { title: yLabel },
      legend: { position: 'none' },
      lineWidth: 2,
      pointSize: 3
    });

  print(title, chart);
}

/* ---------- NBR για Area1 & Area2 ---------- */

// Area1 — NBR
var NBR_A1_IC = yearlyMeanIC('NBR', gArea1, startYear, endYear);
chartIndexSeries(
  NBR_A1_IC,
  'NBR',
  gArea1,
  'Area1 — NBR (annual mean)',
  'NBR',
  30
);

// Area2 — NBR
var NBR_A2_IC = yearlyMeanIC('NBR', gArea2, startYear, endYear);
chartIndexSeries(
  NBR_A2_IC,
  'NBR',
  gArea2,
  'Area2 — NBR (annual mean)',
  'NBR',
  30
);

/* ---------- NDVI for Area3 ---------- */

var NDVI_A3_IC = yearlyMeanIC('NDVI', gArea3, startYear, endYear);
chartIndexSeries(
  NDVI_A3_IC,
  'NDVI',
  gArea3,
  'Area3 — NDVI (annual mean)',
  'NDVI',
  30
);


/***** (SPEI — YEARLY from SPEIbase v2.10) *****/

var speibase12 = ee.ImageCollection('CSIC/SPEI/2_10')
  .filterDate(ee.Date.fromYMD(startYear,1,1), ee.Date.fromYMD(endYear,12,31))
  .select('SPEI_12_month');

function sequenceClamp_(a,b){return ee.List.sequence(ee.Number(a),ee.Number(b));}

function yearlySPEI_mean(ic, geom, y0, y1, label, buffer_m){
  var yearsList = ee.List(sequenceClamp_(y0,y1));
  var geomBuffered = (buffer_m && buffer_m > 0)
    ? ee.Geometry(geom).buffer(buffer_m)
    : ee.Geometry(geom);

  var feats = yearsList.map(function(y){
    y = ee.Number(y);
    var d0 = ee.Date.fromYMD(y,1,1);
    var d1 = d0.advance(1,'year');

    var imYear = ee.Image(ic.filterDate(d0, d1).mean());
    var proj   = imYear.projection();
    var v = imYear.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: ee.Geometry(geomBuffered),
      crs: proj.crs(),
      scale: proj.nominalScale(),
      bestEffort: true,
      maxPixels: 1e9
    }).get('SPEI_12_month');

    return ee.Feature(null, {year: y, spei: v, series: label});
  });

  return ee.FeatureCollection(feats)
           .filter(ee.Filter.notNull(['spei']))
           .sort('year');
}

function yearlySPEI_pixel(ic, geom, y0, y1, label, pixel_buf_m){
  var yearsList = ee.List.sequence(y0, y1);
  var sampleGeom = ee.Geometry(geom).centroid(1).buffer(pixel_buf_m);

  var feats = yearsList.map(function(y){
    y = ee.Number(y);
    var d0 = ee.Date.fromYMD(y,1,1);
    var d1 = d0.advance(1,'year');

    var imYear = ee.Image(ic.filterDate(d0, d1).mean());
    var proj   = imYear.projection();
    var v = imYear.reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: sampleGeom,
      crs: proj.crs(),
      scale: proj.nominalScale(),
      bestEffort: true,
      maxPixels: 1e9
    }).get('SPEI_12_month');

    return ee.Feature(null, {year: y, spei: v, series: label});
  });

  return ee.FeatureCollection(feats)
           .filter(ee.Filter.notNull(['spei']))
           .sort('year');
}

var MEAN_BUFFER_M  = 30000;
var PIXEL_BUFFER_M = 30000;

// Area1
var yA1_mean = yearlySPEI_mean (speibase12, gArea1, startYear, endYear, 'Area1 • mean',  MEAN_BUFFER_M);
var yA1_pix  = yearlySPEI_pixel(speibase12, gArea1, startYear, endYear, 'Area1 • pixel', PIXEL_BUFFER_M);

// Area2
var yA2_mean = yearlySPEI_mean (speibase12, gArea2, startYear, endYear, 'Area2 • mean',  MEAN_BUFFER_M);
var yA2_pix  = yearlySPEI_pixel(speibase12, gArea2, startYear, endYear, 'Area2 • pixel', PIXEL_BUFFER_M);

// Area3
var yA3_mean = yearlySPEI_mean (speibase12, gArea3, startYear, endYear, 'Area3 • mean',  MEAN_BUFFER_M);
var yA3_pix  = yearlySPEI_pixel(speibase12, gArea3, startYear, endYear, 'Area3 • pixel', PIXEL_BUFFER_M);

print('Area1 — Yearly SPEI (mean)',  yA1_mean);
print('Area1 — Yearly SPEI (pixel)', yA1_pix);
print('Area2 — Yearly SPEI (mean)',  yA2_mean);
print('Area2 — Yearly SPEI (pixel)', yA2_pix);
print('Area3 — Yearly SPEI (mean)',  yA3_mean);
print('Area3 — Yearly SPEI (pixel)', yA3_pix);

print(ui.Chart.feature.byFeature(yA1_mean, 'year', ['spei'])
  .setChartType('LineChart')
  .setOptions({title:'Area1 — Yearly SPEI (12m)', hAxis:{title:'Year'}, vAxis:{title:'SPEI'}, lineWidth:2, pointSize:3}));

print(ui.Chart.feature.byFeature(yA2_mean, 'year', ['spei'])
  .setChartType('LineChart')
  .setOptions({title:'Area2 — Yearly SPEI (12m)', hAxis:{title:'Year'}, vAxis:{title:'SPEI'}, lineWidth:2, pointSize:3}));

print(ui.Chart.feature.byFeature(yA3_mean, 'year', ['spei'])
  .setChartType('LineChart')
  .setOptions({title:'Area3 — Yearly SPEI (12m)', hAxis:{title:'Year'}, vAxis:{title:'SPEI'}, lineWidth:2, pointSize:3}));


/***** (4a–4c) Topo + TPI/curv + SOLAR — Copernicus DEM (fixed projection) *****/

/***** Copernicus DEM — single source of truth (NO SRTM) *****/
var DEM_COP_IC = ee.ImageCollection('COPERNICUS/DEM/GLO30').select('DEM');
var DEM_COP    = DEM_COP_IC.mosaic();   // <- single image

// small extent (no global)
var DEM_REGION = ee.Geometry(WE_forestsGeom).buffer(15000).bounds(1);

// Cyprus: UTM 36N ~ EPSG:32636
var WE_DEM_PROJ = ee.Projection('EPSG:32636').atScale(30);

// DEM clipped + fixed projection
var WE_dem = DEM_COP
  .clip(DEM_REGION)
  .setDefaultProjection(WE_DEM_PROJ);

var WE_PI = Math.PI;

// ====== GLOBAL TOPO VARS (for RELAX section) ======
var WE_elev      = WE_dem.rename('elev');
var WE_terr      = ee.Terrain.products(WE_elev);

var WE_slopeDeg  = WE_terr.select('slope');
var WE_aspectDeg = WE_terr.select('aspect');
var WE_aspectRad = WE_aspectDeg.multiply(WE_PI/180);

// (optional)
var WE_TPI300  = WE_elev.subtract(
  WE_elev.focal_mean({kernel: ee.Kernel.square({radius: 300, units:'meters', normalize:true})})
).rename('TPI_300');

var WE_TPI1000 = WE_elev.subtract(
  WE_elev.focal_mean({kernel: ee.Kernel.square({radius: 1000, units:'meters', normalize:true})})
).rename('TPI_1000');

var WE_curv8 = WE_elev.convolve(ee.Kernel.laplacian8()).clamp(-100, 100).rename('curv8');

/* atan(z) for ee.Number μέσω ee.Image.atan + reduceRegion */
function atanNumber_(zNum, geomForEval) {
  var zImg = ee.Image.constant(zNum);
  var angImg = zImg.atan(); // radians
  var pt = ee.Geometry(geomForEval).centroid(1);
  return ee.Number(
    angImg.reduceRegion({
      reducer: ee.Reducer.first(),
      geometry: pt,
      scale: 1000,
      bestEffort: true,
      maxPixels: 1e6
    }).get('constant')
  );
}

/* atan2(y, x) for ee.Number */
function atan2Number_(yNum, xNum, geomForEval) {
  var y = ee.Number(yNum);
  var x = ee.Number(xNum);
  var PI = Math.PI;

  var xIsZero = x.eq(0);
  var yPos    = y.gt(0);
  var yNeg    = y.lt(0);

  var atan_y_over_x = atanNumber_(y.divide(x), geomForEval);

  var angle = ee.Algorithms.If(
    x.gt(0),
    atan_y_over_x,
    ee.Algorithms.If(
      x.lt(0).and(y.gte(0)),
      ee.Number(atan_y_over_x).add(PI),
      ee.Algorithms.If(
        x.lt(0).and(y.lt(0)),
        ee.Number(atan_y_over_x).subtract(PI),
        ee.Algorithms.If(
          xIsZero.and(yPos),
          PI/2,
          ee.Algorithms.If(
            xIsZero.and(yNeg),
            -PI/2,
            0
          )
        )
      )
    )
  );

  return ee.Number(angle);
}

/***** (4a) elev/slope/aspect *****/
function topoSummary(geom, label){
  var geomB = ee.Geometry(geom).buffer(30); // 

var dem = WE_dem.clip(geomB).rename('elev');
  var terrain = ee.Terrain.products(dem);

  var slope   = terrain.select('slope');        // degrees
  var aspectD = terrain.select('aspect');       // degrees 0–360
  var aspectR = aspectD.multiply(Math.PI/180);

  var sinA = aspectR.sin().rename('aspect_sin');
  var cosA = aspectR.cos().rename('aspect_cos');

  var img = dem.addBands([slope, aspectD, sinA, cosA]);

  var stats = img.reduceRegion({
    reducer: ee.Reducer.mean()
      .combine({reducer2: ee.Reducer.stdDev(), sharedInputs:true})
      .combine({reducer2: ee.Reducer.minMax(), sharedInputs:true}),
    geometry: geomB,
    scale: 30,
    bestEffort: true,
    maxPixels: 2e9,
    tileScale: 4
  });

  var d = ee.Dictionary(stats);

  // null-safe circular aspect
  var sinMean = d.get('aspect_sin_mean');
  var cosMean = d.get('aspect_cos_mean');

  var hasVals = ee.Algorithms.If(
    ee.Algorithms.IsEqual(sinMean, null),
    0,
    ee.Algorithms.If(ee.Algorithms.IsEqual(cosMean, null), 0, 1)
  );

  var circMeanRad = ee.Algorithms.If(
    ee.Number(hasVals).eq(1),
    atan2Number_(ee.Number(sinMean), ee.Number(cosMean), geomB),
    null
  );

  var circMeanDeg = ee.Algorithms.If(
    ee.Number(hasVals).eq(1),
    ee.Number(circMeanRad).multiply(180/Math.PI).mod(360),
    null
  );

  var R = ee.Algorithms.If(
    ee.Number(hasVals).eq(1),
    ee.Number(sinMean).pow(2).add(ee.Number(cosMean).pow(2)).sqrt(),
    null
  );

  var circStdDeg = ee.Algorithms.If(
    ee.Number(hasVals).eq(1),
    ee.Number(2).multiply(ee.Number(1).subtract(ee.Number(R))).sqrt().multiply(180/Math.PI),
    null
  );

  return ee.FeatureCollection([ee.Feature(null, {
    area: label,

    elev_mean:   d.get('elev_mean'),
    elev_stdDev: d.get('elev_stdDev'),
    elev_min:    d.get('elev_min'),
    elev_max:    d.get('elev_max'),

    slope_mean_deg:   d.get('slope_mean'),
    slope_stdDev_deg: d.get('slope_stdDev'),
    slope_min_deg:    d.get('slope_min'),
    slope_max_deg:    d.get('slope_max'),

    aspect_mean_deg_raw:   d.get('aspect_mean'),
    aspect_stdDev_deg_raw: d.get('aspect_stdDev'),
    aspect_min_deg:        d.get('aspect_min'),
    aspect_max_deg:        d.get('aspect_max'),

    aspect_circ_mean_deg:  circMeanDeg,
    aspect_resultant_R:    R,
    aspect_circ_std_deg:   circStdDeg,

    northness_mean: ee.Algorithms.If(ee.Number(hasVals).eq(1), cosMean, null),
    eastness_mean:  ee.Algorithms.If(ee.Number(hasVals).eq(1), sinMean, null)
  })]);
}

var topoArea1 = topoSummary(gArea1, 'Area1');
var topoArea2 = topoSummary(gArea2, 'Area2');
var topoArea3 = topoSummary(gArea3, 'Area3');

print('Τοπογραφία — Area1 (elev/slope/aspect):', topoArea1);
print('Τοπογραφία — Area2 (elev/slope/aspect):', topoArea2);
print('Τοπογραφία — Area3 (elev/slope/aspect):', topoArea3);


/***** (4b) TPI_300 / TPI_1000 / curv8 *****/
function tpiCurvSummary(geom, label){
  var geomB = ee.Geometry(geom).buffer(30);
var elev  = WE_dem.clip(geomB).rename('elev');

  var k300  = ee.Kernel.square({radius: 300,  units: 'meters', normalize: true});
  var k1000 = ee.Kernel.square({radius: 1000, units: 'meters', normalize: true});

  var TPI_300  = elev.subtract(elev.focal_mean({kernel: k300})).rename('TPI_300');
  var TPI_1000 = elev.subtract(elev.focal_mean({kernel: k1000})).rename('TPI_1000');
  var curv8    = elev.convolve(ee.Kernel.laplacian8()).clamp(-100, 100).rename('curv8');

  var stack = TPI_300.addBands([TPI_1000, curv8]);

  var stats = stack.reduceRegion({
    reducer: ee.Reducer.mean()
      .combine({reducer2: ee.Reducer.stdDev(), sharedInputs: true})
      .combine({reducer2: ee.Reducer.minMax(), sharedInputs: true}),
    geometry: geomB,
    scale: 30,
    bestEffort: true,
    maxPixels: 2e9,
    tileScale: 4
  });

  var d = ee.Dictionary(stats);

  return ee.FeatureCollection([ee.Feature(null, {
    area: label,
    TPI300_mean: d.get('TPI_300_mean'),   TPI300_stdDev: d.get('TPI_300_stdDev'),
    TPI300_min:  d.get('TPI_300_min'),    TPI300_max:    d.get('TPI_300_max'),
    TPI1000_mean:d.get('TPI_1000_mean'),  TPI1000_stdDev:d.get('TPI_1000_stdDev'),
    TPI1000_min: d.get('TPI_1000_min'),   TPI1000_max:   d.get('TPI_1000_max'),
    curv8_mean:  d.get('curv8_mean'),     curv8_stdDev:  d.get('curv8_stdDev'),
    curv8_min:   d.get('curv8_min'),      curv8_max:     d.get('curv8_max')
  })]);
}

var tpiArea1 = tpiCurvSummary(gArea1, 'Area1');
var tpiArea2 = tpiCurvSummary(gArea2, 'Area2');
var tpiArea3 = tpiCurvSummary(gArea3, 'Area3');

print('TPI/curv8 — Area1:', tpiArea1);
print('TPI/curv8 — Area2:', tpiArea2);
print('TPI/curv8 — Area3:', tpiArea3);


/***** (4c) SOLAR_IDX (APRSEP) *****/
var SOLAR_SEASON = 'APRSEP';
var YEARS_SOLAR  = [1992, 2020];

function _monthsForSeason(mode){
  if (mode === 'DRY')    return [6,7,8,9];
  if (mode === 'APRSEP') return [4,5,6,7,8,9];
  return [1,2,3,4,5,6,7,8,9,10,11,12];
}
var _PI = Math.PI;

var _latR  = ee.Image.pixelLonLat().select('latitude')
  .clip(DEM_REGION)
  .reproject(WE_DEM_PROJ)
  .multiply(_PI/180);


// slope/aspect from DEM_COP
var _terr  = ee.Terrain.products(WE_dem);
var _slope = _terr.select('slope').multiply(_PI/180);
var _aspr  = _terr.select('aspect').multiply(_PI/180);

function solarIndexMonth(year, month){
  var date  = ee.Date.fromYMD(year, month, 15);
  var doy   = date.getRelative('day','year');

  var deltaNum = ee.Number(23.45).multiply(_PI/180)
      .multiply(ee.Number(2*_PI).multiply(doy.add(284)).divide(365).sin());

  var sinδ = ee.Image.constant(deltaNum.sin());
  var cosδ = ee.Image.constant(deltaNum.cos());

  var sinφ = _latR.sin(),      cosφ = _latR.cos();
  var sinβ = _slope.sin(),     cosβ = _slope.cos();
  var cosα = _aspr.cos();

  var cosInc = sinφ.multiply(sinδ).multiply(cosβ)
    .subtract(sinδ.multiply(cosφ).multiply(sinβ).multiply(cosα))
    .add(cosδ.multiply(cosφ).multiply(cosβ))
    .add(cosδ.multiply(sinφ).multiply(sinβ).multiply(cosα));

  return cosInc.max(0).min(1).rename('SOLAR_IDX')
    .set('system:time_start', date.millis())
    .set('year', year).set('month', month);
}

function solarMonthlyIC(y0, y1){
  var years  = ee.List.sequence(y0, y1);
  var months = ee.List.sequence(1,12);
  var imgs = years.map(function(y){
    return months.map(function(m){ return solarIndexMonth(ee.Number(y), ee.Number(m)); });
  }).flatten();
  return ee.ImageCollection(imgs);
}

function solarSeasonForYear(icMonthly, year, seasonMode){
  var months = _monthsForSeason(seasonMode);
  var yearIC = icMonthly.filter(ee.Filter.eq('year', year));
  var seasonIC = yearIC.filter(ee.Filter.inList('month', months));
  return seasonIC.mean().rename('SOLAR_IDX_'+seasonMode).set('year', year);
}

function solarClimatology(icMonthly, y0, y1, seasonMode){
  var years = ee.List.sequence(y0, y1);
  var seasonImgs = years.map(function(y){
    return solarSeasonForYear(icMonthly, ee.Number(y), seasonMode);
  });
  return ee.ImageCollection(seasonImgs).mean()
           .rename('SOLAR_IDX_'+seasonMode+'_clim')
           .set('period', y0 + '-' + y1);
}

function solarStats(img, geom, label){
  var band = ee.String(img.bandNames().get(0));
  var geomB = ee.Geometry(geom).buffer(30);

  var d = ee.Dictionary(img.reduceRegion({
    reducer: ee.Reducer.mean()
      .combine({reducer2: ee.Reducer.stdDev(), sharedInputs:true})
      .combine({reducer2: ee.Reducer.minMax(), sharedInputs:true}),
    geometry: geomB,
    scale: 30,
    bestEffort: true,
    tileScale: 4,
    maxPixels: 2e9
  }));

  return ee.FeatureCollection([ee.Feature(null, {
    area: label, season: SOLAR_SEASON, years: YEARS_SOLAR[0] + '-' + YEARS_SOLAR[1],
    mean:   d.get(band.cat('_mean')),
    stdDev: d.get(band.cat('_stdDev')),
    min:    d.get(band.cat('_min')),
    max:    d.get(band.cat('_max'))
  })]);
}

var solarMonthly = solarMonthlyIC(YEARS_SOLAR[0], YEARS_SOLAR[1]);
var solarClim    = solarClimatology(solarMonthly, YEARS_SOLAR[0], YEARS_SOLAR[1], SOLAR_SEASON);

var solArea1 = solarStats(solarClim, gArea1, 'Area1');
var solArea2 = solarStats(solarClim, gArea2, 'Area2');
var solArea3 = solarStats(solarClim, gArea3, 'Area3');

print('SOLAR_IDX — Area1 ('+SOLAR_SEASON+'):', solArea1);
print('SOLAR_IDX — Area2 ('+SOLAR_SEASON+'):', solArea2);
print('SOLAR_IDX — Area3 ('+SOLAR_SEASON+'):', solArea3);



/***** (GEOLOGY — minimal, using field NAME) Dominant geology by area ****/

var GEOL_FC = ee.FeatureCollection('users/filipposeliades31/GEOLOGY');
var GEOL_FIELD = 'NAME';

var cats = ee.List(GEOL_FC.aggregate_array(GEOL_FIELD)).distinct().sort();

var GEOL_CODED = GEOL_FC.map(function (f) {
  var name = f.get(GEOL_FIELD);
  var code = cats.indexOf(name);
  return f.set('code', code);
});

var geolCodeImg = GEOL_CODED.reduceToImage({
  properties: ['code'],
  reducer: ee.Reducer.first()
}).rename('geol_code');

function geologyTop(geom, label, scale) {
  scale = scale || 30;
  var stack = ee.Image.pixelArea().rename('area').addBands(geolCodeImg).clip(geom);

  var grouped = stack.reduceRegion({
    reducer: ee.Reducer.sum().group({ groupField: 1, groupName: 'code' }),
    geometry: geom, scale: scale, bestEffort: true, tileScale: 4, maxPixels: 2e9
  });

  var groups = ee.List(
    ee.Algorithms.If(
      ee.Dictionary(grouped).contains('groups'),
      ee.Dictionary(grouped).get('groups'),
      ee.List([])
    )
  );

  var fc = ee.FeatureCollection(groups.map(function (g) {
    g = ee.Dictionary(g);
    var code = ee.Number(g.get('code')).toInt();
    var a_m2 = ee.Number(g.get('sum'));
    var name = ee.String(cats.get(code));
    return ee.Feature(null, {
      area: label,
      code: code,
      name: name,
      area_km2: a_m2.divide(1e6)
    });
  })).sort('area_km2', false);

  return fc;
}

var topArea1 = geologyTop(gArea1, 'Area1', 30);
var topArea2 = geologyTop(gArea2, 'Area2', 30);
var topArea3 = geologyTop(gArea3, 'Area3', 30);

print('Geology — Area1 (Top-3):', topArea1.limit(3));
print('Geology — Area2 (Top-3):', topArea2.limit(3));
print('Geology — Area3 (Top-3):', topArea3.limit(3));

/***** ============================================================
WHERE ELSE
============================================================ *****/


var WE_SEASON_MODE = 'APRSEP';
var WE_EVENT_YEARS = { Area1: 2016, Area2: 2008, Area3: 2008 };

// Elevation toggle
var USE_ELEVATION = false;

// SPEI threshold
var WE_SPEI_THR_ALL = -1.4;

// =====================
// RELAX tuning (μόνο RELAX)
// =====================
// Event DROP% relaxation:
// 1.00 = none, 0.90 slight, 0.80 medium, 0.65 strong
// (if the threshold is negative, factor < 1 => closer to 0 => more relaxed)
var RELAX_DROP_FACTOR = 0.60;

// Do I want to require persistence +1 year in RELAX?
var REQUIRE_PERSIST_POST1 = true;

// Persistence relaxation in RELAX, but on DELTA (target - baseline):
// 1.00 = none, 0.90 slight, 0.80 medium, 0.65 strong
var RELAX_PERSIST_FACTOR = 0.50;

// OPTIONAL: minimum “negative trend” strength at post+1 as a fraction of the EVENT delta.
// - Set 0.00 if you want ONLY: (DELTA_post1 < 0)
// - 0.30 = low strength, 0.50 = medium, 0.70 = stricter
var PERSIST_DELTA_FRAC = 0.00;


function WE_monthsForSeason(mode){
return ee.List(
mode === 'DRY'    ? [6,7,8,9] :
mode === 'APRSEP' ? [4,5,6,7,8,9] :
[1,2,3,4,5,6,7,8,9,10,11,12]
);
}

function WE_aspectWindowDeg(centerDeg, halfWidthDeg){
var mu  = ee.Image.constant(ee.Number(centerDeg).multiply(Math.PI/180));
var wid = ee.Image.constant(ee.Number(halfWidthDeg).max(1).multiply(Math.PI/180));
var cosw = wid.cos();
var cosDiff = WE_aspectRad.subtract(mu).cos();
return cosDiff.gte(cosw);
}

// Geology masks/helpers 
function geologyIncludeMask(codesList){
var init = ee.Image(0).rename('geol_inc');
var mask = ee.List(codesList).iterate(function(code, acc){
acc = ee.Image(acc);
var m = geolCodeImg.eq(ee.Number(code)).rename('m_code');
return acc.or(m);
}, init);
return ee.Image(mask);
}
function geologyTopCodes(areaLabel, topN){
var fc = (areaLabel === 'Area1') ? topArea1
: (areaLabel === 'Area2') ? topArea2
: topArea3;
var fcTop = fc.limit(topN);
return ee.List(fcTop.aggregate_array('code'));
}

// Elevation mask from stats 
var ELEV_STD_MULT = 1.0;
function elevationMaskFromStats(areaLabel){
var d = ee.Dictionary(
areaLabel === 'Area1' ? topoArea1.first().toDictionary()
: areaLabel === 'Area2' ? topoArea2.first().toDictionary()
: topoArea3.first().toDictionary()
);
var mu  = ee.Number(d.get('elev_mean'));
var sd  = ee.Number(d.get('elev_stdDev'));
var lo  = mu.subtract(sd.multiply(ELEV_STD_MULT)).max(0);
var hi  = mu.add(sd.multiply(ELEV_STD_MULT));
return {
mask: WE_elev.gte(lo).and(WE_elev.lte(hi)),
lo: lo, hi: hi
};
}

// Solar index (season mean) per year
var WE_latRad = ee.Image.pixelLonLat().select('latitude')
  .clip(DEM_REGION)
  .reproject(WE_DEM_PROJ)
  .multiply(WE_PI/180);
var WE_slopeRad = WE_slopeDeg.multiply(WE_PI/180);

function WE_solarIndexMonth(y, m){
var date  = ee.Date.fromYMD(y, m, 15);
var doy   = date.getRelative('day','year');
var delta = ee.Number(23.45).multiply(WE_PI/180)
.multiply(ee.Number(2*WE_PI).multiply(doy.add(284)).divide(365).sin());
var sinδ = ee.Image.constant(delta.sin()), cosδ = ee.Image.constant(delta.cos());
var sinφ = WE_latRad.sin(),  cosφ = WE_latRad.cos();
var sinβ = WE_slopeRad.sin(), cosβ = WE_slopeRad.cos();
var cosα = WE_aspectRad.cos();
var cosInc = sinφ.multiply(sinδ).multiply(cosβ)
.subtract(sinδ.multiply(cosφ).multiply(sinβ).multiply(cosα))
.add(cosδ.multiply(cosφ).multiply(cosβ))
.add(cosδ.multiply(sinφ).multiply(sinβ).multiply(cosα));
return cosInc.max(0).min(1).rename('SOLAR_IDX')
.set('year', y).set('month', m).set('system:time_start', date.millis());
}
function WE_solarSeasonYear(y, season){
var mm = WE_monthsForSeason(season);
var ic = ee.ImageCollection(mm.map(function(m){ return WE_solarIndexMonth(y, ee.Number(m)); }));
return ic.mean().rename('SOLAR_IDX');
}

// SPEI
var WE_speibase12 = ee.ImageCollection('CSIC/SPEI/2_10').select('SPEI_12_month');
function WE_speiSeasonYear(y, season){
var mm  = WE_monthsForSeason(season);
var d0  = ee.Date.fromYMD(y,1,1), d1 = d0.advance(1,'year');
var icY = WE_speibase12.filterDate(d0, d1)
.map(function(im){ return im.set('m', ee.Date(im.get('system:time_start')).get('month')); })
.filter(ee.Filter.inList('m', mm));
return icY.mean().rename('SPEI_12m');
}

// ============================================================
// TRUTH masks = AOI ∩ forest
// ============================================================
var A1_mortMask = ee.Image(1).clip(gArea1).updateMask(WE_forestMask30).rename('A1_mort');
var A2_mortMask = ee.Image(1).clip(gArea2).updateMask(WE_forestMask30).rename('A2_mort');
var A3_mortMask = ee.Image(1).clip(gArea3).updateMask(WE_forestMask30).rename('A3_mort');

// ============================================================
// METRICS vs pre3: base = mean(y-3..y-1)
// Returns 2 bands:
//  - <band>_drop  (% change vs base)
//  - <band>_delta (target - base)  --> persistence criterion
// NOTE: fixed band names (no ee.Number inside the name) => no rename errors
// ============================================================

function WE_metricsVsPre3(band, roi, eventYear, targetYear){
var base = ee.ImageCollection(
ee.List.sequence(eventYear-3, eventYear-1).map(function(y){
return WE_yearMeanImage(band, y, roi);
})
).mean();

var tgt  = WE_yearMeanImage(band, targetYear, roi);

var delta = tgt.subtract(base).rename(band + '_delta');
var drop  = delta.divide(base).multiply(100).rename(band + '_drop');

return drop.addBands(delta);
}

// ----- EVENT (target=eventYear) -----
var MET_NDVI_EVT_A1 = WE_metricsVsPre3('NDVI', WE_forestsGeom, WE_EVENT_YEARS.Area1, WE_EVENT_YEARS.Area1);
var MET_NBR_EVT_A1  = WE_metricsVsPre3('NBR',  WE_forestsGeom, WE_EVENT_YEARS.Area1, WE_EVENT_YEARS.Area1);

var MET_NDVI_EVT_A2 = WE_metricsVsPre3('NDVI', WE_forestsGeom, WE_EVENT_YEARS.Area2, WE_EVENT_YEARS.Area2);
var MET_NBR_EVT_A2  = WE_metricsVsPre3('NBR',  WE_forestsGeom, WE_EVENT_YEARS.Area2, WE_EVENT_YEARS.Area2);

var MET_NDVI_EVT_A3 = WE_metricsVsPre3('NDVI', WE_forestsGeom, WE_EVENT_YEARS.Area3, WE_EVENT_YEARS.Area3);
var MET_NBR_EVT_A3  = WE_metricsVsPre3('NBR',  WE_forestsGeom, WE_EVENT_YEARS.Area3, WE_EVENT_YEARS.Area3);

// DROP% images (event)
var DROP_NDVI_A1 = MET_NDVI_EVT_A1.select('NDVI_drop').rename('NDVI_drop_evt_A1');
var DROP_NBR_A1  = MET_NBR_EVT_A1.select('NBR_drop').rename('NBR_drop_evt_A1');

var DROP_NDVI_A2 = MET_NDVI_EVT_A2.select('NDVI_drop').rename('NDVI_drop_evt_A2');
var DROP_NBR_A2  = MET_NBR_EVT_A2.select('NBR_drop').rename('NBR_drop_evt_A2');

var DROP_NDVI_A3 = MET_NDVI_EVT_A3.select('NDVI_drop').rename('NDVI_drop_evt_A3');
var DROP_NBR_A3  = MET_NBR_EVT_A3.select('NBR_drop').rename('NBR_drop_evt_A3');

// EVENT delta images (για optional δέσιμο με PERSIST_DELTA_FRAC)
var DELTA_NDVI_EVT_A1 = MET_NDVI_EVT_A1.select('NDVI_delta').rename('NDVI_delta_evt_A1');
var DELTA_NBR_EVT_A1  = MET_NBR_EVT_A1.select('NBR_delta').rename('NBR_delta_evt_A1');

var DELTA_NDVI_EVT_A2 = MET_NDVI_EVT_A2.select('NDVI_delta').rename('NDVI_delta_evt_A2');
var DELTA_NBR_EVT_A2  = MET_NBR_EVT_A2.select('NBR_delta').rename('NBR_delta_evt_A2');

var DELTA_NDVI_EVT_A3 = MET_NDVI_EVT_A3.select('NDVI_delta').rename('NDVI_delta_evt_A3');
var DELTA_NBR_EVT_A3  = MET_NBR_EVT_A3.select('NBR_delta').rename('NBR_delta_evt_A3');

// ----- POST+1 (target=eventYear+1) -----
var Y1_A1 = ee.Number(WE_EVENT_YEARS.Area1).add(1);
var Y1_A2 = ee.Number(WE_EVENT_YEARS.Area2).add(1);
var Y1_A3 = ee.Number(WE_EVENT_YEARS.Area3).add(1);

var MET_NDVI_P1_A1 = WE_metricsVsPre3('NDVI', WE_forestsGeom, WE_EVENT_YEARS.Area1, Y1_A1);
var MET_NBR_P1_A1  = WE_metricsVsPre3('NBR',  WE_forestsGeom, WE_EVENT_YEARS.Area1, Y1_A1);

var MET_NDVI_P1_A2 = WE_metricsVsPre3('NDVI', WE_forestsGeom, WE_EVENT_YEARS.Area2, Y1_A2);
var MET_NBR_P1_A2  = WE_metricsVsPre3('NBR',  WE_forestsGeom, WE_EVENT_YEARS.Area2, Y1_A2);

var MET_NDVI_P1_A3 = WE_metricsVsPre3('NDVI', WE_forestsGeom, WE_EVENT_YEARS.Area3, Y1_A3);
var MET_NBR_P1_A3  = WE_metricsVsPre3('NBR',  WE_forestsGeom, WE_EVENT_YEARS.Area3, Y1_A3);

// PERSIST delta images (post+1): target - baseline_pre3
var DELTA_NDVI_P1_A1 = MET_NDVI_P1_A1.select('NDVI_delta').rename('NDVI_delta_post1_A1');
var DELTA_NBR_P1_A1  = MET_NBR_P1_A1.select('NBR_delta').rename('NBR_delta_post1_A1');

var DELTA_NDVI_P1_A2 = MET_NDVI_P1_A2.select('NDVI_delta').rename('NDVI_delta_post1_A2');
var DELTA_NBR_P1_A2  = MET_NBR_P1_A2.select('NBR_delta').rename('NBR_delta_post1_A2');

var DELTA_NDVI_P1_A3 = MET_NDVI_P1_A3.select('NDVI_delta').rename('NDVI_delta_post1_A3');
var DELTA_NBR_P1_A3  = MET_NBR_P1_A3.select('NBR_delta').rename('NBR_delta_post1_A3');

// ============================================================
// AUTO thresholds από TRUTH AOI
// - Event thresholds σε %DROP 
// - Persistence thresholds in DELTA (for a “negative trend” regardless of the % sign)
// ============================================================
var WE_MIN_DROP_NBR;
var WE_MIN_DROP_NDVI;

var WE_EVT_DELTA_NBR;
var WE_EVT_DELTA_NDVI;

var WE_P1_DELTA_NBR;
var WE_P1_DELTA_NDVI;

function meanInTruth(singleBandImg, roiGeom, mortMask, label){
var geom = ee.Geometry(roiGeom).intersection(WE_forestsGeom, 1);

var v = ee.Image(singleBandImg)
.updateMask(WE_forestMask30)
.updateMask(mortMask)
.reduceRegion({
reducer: ee.Reducer.mean(),
geometry: geom,
scale: 30,
bestEffort: true,
tileScale: 4,
maxPixels: 2e9
}).values().get(0);

v = ee.Number(v);
//print('--- ' + label + ' mean inside TRUTH (AOI∩forest) ---', v);
return v;
}

// --- Event %drop thresholds from truth ---
var MIN_DROP_NBR_A1  = meanInTruth(DROP_NBR_A1,  gArea1, A1_mortMask, 'Area1 [NBR %DROP evt]');
var MIN_DROP_NBR_A2  = meanInTruth(DROP_NBR_A2,  gArea2, A2_mortMask, 'Area2 [NBR %DROP evt]');
var MIN_DROP_NDVI_A3 = meanInTruth(DROP_NDVI_A3, gArea3, A3_mortMask, 'Area3 [NDVI %DROP evt]');

WE_MIN_DROP_NBR  = ee.Number(MIN_DROP_NBR_A1).add(MIN_DROP_NBR_A2).divide(2);
WE_MIN_DROP_NDVI = ee.Number(MIN_DROP_NDVI_A3);

// Leave-one-AOI-out//
//WE_MIN_DROP_NBR  = ee.Number(MIN_DROP_NBR_A2); 
//WE_MIN_DROP_NDVI = ee.Number(MIN_DROP_NDVI_A3);

//WE_MIN_DROP_NBR  = ee.Number(MIN_DROP_NBR_A1); 
//WE_MIN_DROP_NDVI = ee.Number(MIN_DROP_NDVI_A3); 

// --- Event delta thresholds from truth 
var EVT_DELTA_NBR_A1  = meanInTruth(DELTA_NBR_EVT_A1,  gArea1, A1_mortMask, 'Area1 [NBR DELTA evt]');
var EVT_DELTA_NBR_A2  = meanInTruth(DELTA_NBR_EVT_A2,  gArea2, A2_mortMask, 'Area2 [NBR DELTA evt]');
var EVT_DELTA_NDVI_A3 = meanInTruth(DELTA_NDVI_EVT_A3, gArea3, A3_mortMask, 'Area3 [NDVI DELTA evt]');

WE_EVT_DELTA_NBR  = ee.Number(EVT_DELTA_NBR_A1).add(EVT_DELTA_NBR_A2).divide(2);
WE_EVT_DELTA_NDVI = ee.Number(EVT_DELTA_NDVI_A3);

// --- Post+1 persistence delta thresholds from truth ---
var P1_DELTA_NBR_A1  = meanInTruth(DELTA_NBR_P1_A1,  gArea1, A1_mortMask, 'Area1 [NBR DELTA post+1]');
var P1_DELTA_NBR_A2  = meanInTruth(DELTA_NBR_P1_A2,  gArea2, A2_mortMask, 'Area2 [NBR DELTA post+1]');
var P1_DELTA_NDVI_A3 = meanInTruth(DELTA_NDVI_P1_A3, gArea3, A3_mortMask, 'Area3 [NDVI DELTA post+1]');

WE_P1_DELTA_NBR  = ee.Number(P1_DELTA_NBR_A1).add(P1_DELTA_NBR_A2).divide(2);
WE_P1_DELTA_NDVI = ee.Number(P1_DELTA_NDVI_A3);

// console summary
//print('================ FINAL AUTO THRESHOLDS (from AOI truth) ================');
//print('WE_MIN_DROP_NBR   (event %DROP):', WE_MIN_DROP_NBR);
//print('WE_MIN_DROP_NDVI  (event %DROP):', WE_MIN_DROP_NDVI);
//print('WE_EVT_DELTA_NBR  (event DELTA):', WE_EVT_DELTA_NBR);
//print('WE_EVT_DELTA_NDVI (event DELTA):', WE_EVT_DELTA_NDVI);
//print('WE_P1_DELTA_NBR   (post+1 DELTA):', WE_P1_DELTA_NBR);
//print('WE_P1_DELTA_NDVI  (post+1 DELTA):', WE_P1_DELTA_NDVI);
//print('RELAX_DROP_FACTOR:', RELAX_DROP_FACTOR,
//' | REQUIRE_PERSIST_POST1:', REQUIRE_PERSIST_POST1,
//' | RELAX_PERSIST_FACTOR:', RELAX_PERSIST_FACTOR,
//' | PERSIST_DELTA_FRAC:', PERSIST_DELTA_FRAC);

print('===== THRESHOLDS USED (EVENT %DROP) =====');
print('NBR STRICT thr_event_drop_pct =', WE_MIN_DROP_NBR);
print('NBR RELAX  thr_event_drop_pct =', ee.Number(WE_MIN_DROP_NBR).multiply(RELAX_DROP_FACTOR));

print('NDVI STRICT thr_event_drop_pct =', WE_MIN_DROP_NDVI);
print('NDVI RELAX  thr_event_drop_pct =', ee.Number(WE_MIN_DROP_NDVI).multiply(RELAX_DROP_FACTOR));
print('========================================');


print('========================================================================');

// ============================================================
// RELAXED TOPO CONTROLS FOR LIKE MASKS
// - NDVI/NBR = EVENT %DROP + (RELAX only) PERSIST via DELTA(post+1)
// - Persistence = “negative trend” => DELTA_post1 < 0
//   + optional knob: PERSIST_DELTA_FRAC (0 => only <0)
// ============================================================

/***** RELAX PERCENTAGE *****/
var MAJORITY_FACTOR = 0.6; // 

// Topographic toggles
var USE_TOPO     = true;
var USE_ASPECT   = true;
var USE_SLOPE    = true;
var USE_TPI300   = true;
var USE_TPI1000  = true;
var USE_CURV8    = true;

// Other criteria
var USE_SOLAR     = true;
var USE_SPEI      = true;
var USE_GEOLOGY   = true;

// Tele indices by area (as you want it)
var USE_NDVI_like = { Area1: false, Area2: false, Area3: true  };
var USE_NBR_like  = { Area1: true,  Area2: true,  Area3: false };


// Strict/Relax
var STRICT_K    = 1;
var STRICT_RULE = 'AND';
var RELAX_K     = 1.5;
var RELAX_RULE  = 'MAJORITY';

// Quantiles
var Q_STRICT_LOW  = 0.25;
var Q_STRICT_HIGH = 0.75;
var Q_RELAX_LOW   = 0.05;
var Q_RELAX_HIGH  = 0.95;

// Helpers
function _clamp(v, lo, hi){ return ee.Number(v).max(lo).min(hi); }
function _wrapDeg(x){ return ee.Number(x).mod(360); }
function _num(d, k, fb){ d = ee.Dictionary(d); return ee.Number(ee.Algorithms.If(d.contains(k), d.get(k), fb)); }
function _toBinBand(m){ return ee.Image(m).gt(0).unmask(0).rename('b').toByte(); }
function _bandRangeMask(img, lo, hi){ return ee.Image(img).gte(lo).and(ee.Image(img).lte(hi)); }

function _combineMasks(masks, rule){
masks = ee.List(masks);
var coll = ee.ImageCollection(masks.map(function(m){ return _toBinBand(m); }));
var N    = coll.size();

var emptyIsNeutral = ee.Algorithms.If(
N.eq(0),
ee.Image(1).rename('b').toByte(),
(function(){
var sum  = coll.sum();
var need = ee.Number(N).multiply(MAJORITY_FACTOR).ceil(); // 60% majority
return ee.Image(ee.Algorithms.If(
rule === 'AND', sum.eq(N), sum.gte(need)
)).rename('b').toByte();
})()
);

return ee.Image(emptyIsNeutral);
}

function getQuantileRange(img, roiGeom, qLow, qHigh){
  var band = ee.String(img.bandNames().get(0));

  var pLow  = ee.Number(qLow).multiply(100);
  var pHigh = ee.Number(qHigh).multiply(100);

  var qdict = ee.Dictionary(
    img.updateMask(WE_forestMask30).reduceRegion({
      reducer: ee.Reducer.percentile([pLow, pHigh]),
      geometry: roiGeom,
      scale: 90,
      bestEffort: true,
      tileScale: 4,
      maxPixels: 2e9
    })
  );

  var keyLo = band.cat('_p').cat(pLow.format('%.0f'));
  var keyHi = band.cat('_p').cat(pHigh.format('%.0f'));

  var lo = ee.Number(ee.Algorithms.If(
    qdict.contains(keyLo),
    qdict.get(keyLo),
    img.reduceRegion({
      reducer: ee.Reducer.min(),
      geometry: roiGeom,
      scale: 90,
      bestEffort: true,
      tileScale: 4,
      maxPixels: 2e9
    }).get(band)
  ));

  var hi = ee.Number(ee.Algorithms.If(
    qdict.contains(keyHi),
    qdict.get(keyHi),
    img.reduceRegion({
      reducer: ee.Reducer.max(),
      geometry: roiGeom,
      scale: 90,
      bestEffort: true,
      tileScale: 4,
      maxPixels: 2e9
    }).get(band)
  ));

  return { lo: lo, hi: hi };
}


// SPEI timing options
var SPEI_MODE = 'EVENT_MEAN';

function WE_speiHydroYear(y, startMonth) {
var d0 = ee.Date.fromYMD(y-1, startMonth, 1);
var d1 = ee.Date.fromYMD(y,   startMonth, 1);
return WE_speibase12.filterDate(d0, d1).mean().rename('SPEI_12m');
}
function WE_speiSeasonYear_min(y, season){
var mm = WE_monthsForSeason(season);
var d0 = ee.Date.fromYMD(y,1,1), d1 = d0.advance(1,'year');
var icY = WE_speibase12.filterDate(d0, d1)
.map(function(im){ return im.set('m', ee.Date(im.get('system:time_start')).get('month')); })
.filter(ee.Filter.inList('m', mm));
return icY.min().rename('SPEI_12m_min');
}
function WE_speiCalendarYear(y){
var d0 = ee.Date.fromYMD(y, 1, 1);
var d1 = d0.advance(1, 'year');
return WE_speibase12.filterDate(d0, d1).mean().rename('SPEI_12m');
}
function getSpeiMask(areaLabel){
var y = WE_EVENT_YEARS[areaLabel];
var speiImg = ee.Image(1);

speiImg = ee.Image(ee.Algorithms.If(
SPEI_MODE === 'EVENT_MEAN',
WE_speiSeasonYear(y, WE_SEASON_MODE),
speiImg
));
speiImg = ee.Image(ee.Algorithms.If(
SPEI_MODE === 'EVENT_MIN',
WE_speiSeasonYear_min(y, WE_SEASON_MODE),
speiImg
));
speiImg = ee.Image(ee.Algorithms.If(
SPEI_MODE === 'LAG_MEAN',
WE_speiSeasonYear(y-1, WE_SEASON_MODE),
speiImg
));
speiImg = ee.Image(ee.Algorithms.If(
SPEI_MODE === 'HYDRO_MEAN',
WE_speiHydroYear(y, 10),
speiImg
));
speiImg = ee.Image(ee.Algorithms.If(
SPEI_MODE === 'CAL_MEAN',
WE_speiCalendarYear(y),
speiImg
));

return speiImg.lte(WE_SPEI_THR_ALL);
}

// -------------------------
// makeLikeMask (EVENT %DROP + RELAX-only PERSIST via DELTA)
// -------------------------
function makeLikeMask(areaLabel, K, combineRule, modeName){
var isRelax = (modeName === 'RELAX');

var roiGeom = (areaLabel === 'Area1') ? gArea1
: (areaLabel === 'Area2') ? gArea2
: gArea3;

var qLow  = (modeName === 'STRICT') ? Q_STRICT_LOW  : Q_RELAX_LOW;
var qHigh = (modeName === 'STRICT') ? Q_STRICT_HIGH : Q_RELAX_HIGH;

var topoD = ee.Dictionary(
areaLabel==='Area1' ? topoArea1.first().toDictionary()
: areaLabel==='Area2' ? topoArea2.first().toDictionary()
: topoArea3.first().toDictionary()
);
var solD  = ee.Dictionary(
areaLabel==='Area1' ? solArea1.first().toDictionary()
: areaLabel==='Area2' ? solArea2.first().toDictionary()
: solArea3.first().toDictionary()
);

// ----- TOPO masks -----
var topoMasks = [];

if (USE_TOPO && USE_ASPECT){
var aspect_mu = _wrapDeg(_num(topoD,'aspect_circ_mean_deg',180));
var aspect_sd = _num(topoD,'aspect_circ_std_deg',30).multiply(K);
topoMasks.push(WE_aspectWindowDeg(aspect_mu, aspect_sd));
}

if (USE_TOPO && USE_SLOPE){
var slopeQ = getQuantileRange(WE_slopeDeg, roiGeom, qLow, qHigh);
var slope_lo  = _clamp(slopeQ.lo, 0, 90);
var slope_hi  = _clamp(slopeQ.hi, 0, 90);
topoMasks.push(WE_slopeDeg.gte(slope_lo).and(WE_slopeDeg.lte(slope_hi)));
}

if (USE_TOPO && USE_TPI300){
var tpi3Q = getQuantileRange(WE_TPI300, roiGeom, qLow, qHigh);
topoMasks.push(_bandRangeMask(WE_TPI300, tpi3Q.lo, tpi3Q.hi));
}

if (USE_TOPO && USE_TPI1000){
var tpi1Q = getQuantileRange(WE_TPI1000, roiGeom, qLow, qHigh);
topoMasks.push(_bandRangeMask(WE_TPI1000, tpi1Q.lo, tpi1Q.hi));
}

if (USE_TOPO && USE_CURV8){
var curvQ = getQuantileRange(WE_curv8, roiGeom, qLow, qHigh);
topoMasks.push(_bandRangeMask(WE_curv8, curvQ.lo, curvQ.hi));
}

var m_topo = _combineMasks(topoMasks, combineRule);

// ----- SOLAR -----
var solar_mu  = _num(solD,'mean',0.90);
var solar_sd  = _num(solD,'stdDev',0.03);
var solar_min = ee.Number(solar_mu).subtract(solar_sd).max(0).min(1);
var useSolarThisArea = (areaLabel!=='Area2') && USE_SOLAR;

var m_solar = useSolarThisArea
? WE_solarSeasonYear(WE_EVENT_YEARS[areaLabel], WE_SEASON_MODE).gte(solar_min)
: ee.Image(1);

// ----- SPEI -----
var m_spei = USE_SPEI ? getSpeiMask(areaLabel) : ee.Image(1);

// ----- Elevation -----
var elevObj = USE_ELEVATION ? elevationMaskFromStats(areaLabel) : {mask: ee.Image(1)};
var m_elev  = ee.Image(elevObj.mask);

// ----- Geology -----
var geocodes = geologyTopCodes(areaLabel, 3);
var m_geol   = USE_GEOLOGY ? geologyIncludeMask(geocodes) : ee.Image(1);

// =========================================================
// NDVI criterion (Area3 only): EVENT %DROP + (RELAX only) PERSIST DELTA
// =========================================================
var m_ndvi;
if (USE_NDVI_like[areaLabel]) {

var dropEvt = (areaLabel === 'Area1') ? DROP_NDVI_A1
           : (areaLabel === 'Area2') ? DROP_NDVI_A2
                                     : DROP_NDVI_A3;

var deltaEvt = (areaLabel === 'Area1') ? DELTA_NDVI_EVT_A1
            : (areaLabel === 'Area2') ? DELTA_NDVI_EVT_A2
                                      : DELTA_NDVI_EVT_A3;

var deltaP1  = (areaLabel === 'Area1') ? DELTA_NDVI_P1_A1
            : (areaLabel === 'Area2') ? DELTA_NDVI_P1_A2
                                      : DELTA_NDVI_P1_A3;

// Event threshold (%DROP): strict = truth, relax = relaxed
var thrEvt = ee.Number(WE_MIN_DROP_NDVI);
var thrUsed = ee.Number(ee.Algorithms.If(isRelax, thrEvt.multiply(RELAX_DROP_FACTOR), thrEvt));
var m_evt = dropEvt.lte(thrUsed);

// Persistence threshold (DELTA):
// - BASIC: deltaP1 < 0  (negative trend)
// - + optional minimum intensity:
//   deltaP1 <= max( truth_post1 * (relaxFactor), eventDelta * PERSIST_DELTA_FRAC, 0 )
//   (the max pulls the threshold closer to 0 => more relaxed)

var thrTruthP1 = ee.Number(WE_P1_DELTA_NDVI).multiply(isRelax ? RELAX_PERSIST_FACTOR : 1.0);
var thrFrac    = ee.Number(WE_EVT_DELTA_NDVI).multiply(PERSIST_DELTA_FRAC);
var thrP1      = ee.Number(thrTruthP1).max(thrFrac).max(0); // αν PERSIST_DELTA_FRAC=0 => γίνεται 0

var m_p1 = deltaP1.lt(0).and(deltaP1.lte(thrP1));

m_ndvi = ee.Image(ee.Algorithms.If(
  isRelax && REQUIRE_PERSIST_POST1,
  m_evt.and(m_p1),
  m_evt
));


} else {
m_ndvi = ee.Image(1);
}

// =========================================================
// NBR criterion (Area1/2): EVENT %DROP + (RELAX only) PERSIST DELTA
// =========================================================
var m_nbr;
if (USE_NBR_like[areaLabel]) {


var dropEvtN = (areaLabel === 'Area1') ? DROP_NBR_A1
            : (areaLabel === 'Area2') ? DROP_NBR_A2
                                      : DROP_NBR_A3;

var deltaEvtN = (areaLabel === 'Area1') ? DELTA_NBR_EVT_A1
             : (areaLabel === 'Area2') ? DELTA_NBR_EVT_A2
                                       : DELTA_NBR_EVT_A3;

var deltaP1N  = (areaLabel === 'Area1') ? DELTA_NBR_P1_A1
             : (areaLabel === 'Area2') ? DELTA_NBR_P1_A2
                                       : DELTA_NBR_P1_A3;

// Event threshold (%DROP)
var thrEvtN = ee.Number(WE_MIN_DROP_NBR);
var thrUsedN = ee.Number(ee.Algorithms.If(isRelax, thrEvtN.multiply(RELAX_DROP_FACTOR), thrEvtN));
var m_evtN = dropEvtN.lte(thrUsedN);

// Persistence threshold (DELTA)
var thrTruthP1N = ee.Number(WE_P1_DELTA_NBR).multiply(isRelax ? RELAX_PERSIST_FACTOR : 1.0);
var thrFracN    = ee.Number(WE_EVT_DELTA_NBR).multiply(PERSIST_DELTA_FRAC);
var thrP1N      = ee.Number(thrTruthP1N).max(thrFracN).max(0);

var m_p1N = deltaP1N.lt(0).and(deltaP1N.lte(thrP1N));

m_nbr = ee.Image(ee.Algorithms.If(
  isRelax && REQUIRE_PERSIST_POST1,
  m_evtN.and(m_p1N),
  m_evtN
));


} else {
m_nbr = ee.Image(1);
}

return m_topo.and(m_solar).and(m_spei).and(m_ndvi).and(m_nbr).and(m_elev).and(m_geol)
.updateMask(WE_forestMask30)
.rename('b');
}

// STRICT vs RELAX masks
var A1_like_STRICT = makeLikeMask('Area1', STRICT_K, STRICT_RULE, 'STRICT').rename('A1_like_STRICT');
var A2_like_STRICT = makeLikeMask('Area2', STRICT_K, STRICT_RULE, 'STRICT').rename('A2_like_STRICT');
var A3_like_STRICT = makeLikeMask('Area3', STRICT_K, STRICT_RULE, 'STRICT').rename('A3_like_STRICT');

var A1_like_RELAX  = makeLikeMask('Area1', RELAX_K, RELAX_RULE, 'RELAX').rename('A1_like_RELAX');
var A2_like_RELAX  = makeLikeMask('Area2', RELAX_K, RELAX_RULE, 'RELAX').rename('A2_like_RELAX');
var A3_like_RELAX  = makeLikeMask('Area3', RELAX_K, RELAX_RULE, 'RELAX').rename('A3_like_RELAX');

// Layers
Map.addLayer(A1_like_STRICT.selfMask(), {min:0,max:1,palette:['#fb6a4a']}, 'Area1 — STRICT topo (quantiles)', false);
Map.addLayer(A1_like_RELAX.selfMask(),  {min:0,max:1,palette:['#cb181d']}, 'Area1 — RELAXED topo (quantiles)', true);

Map.addLayer(A2_like_STRICT.selfMask(), {min:0,max:1,palette:['#6baed6']}, 'Area2 — STRICT topo (quantiles)', false);
Map.addLayer(A2_like_RELAX.selfMask(),  {min:0,max:1,palette:['#08519c']}, 'Area2 — RELAXED topo (quantiles)', false);

Map.addLayer(A3_like_STRICT.selfMask(), {min:0,max:1,palette:['#74c476']}, 'Area3 — STRICT topo (quantiles)', false);
Map.addLayer(A3_like_RELAX.selfMask(),  {min:0,max:1,palette:['#006d2c']}, 'Area3 — RELAXED topo (quantiles)', false);

// Export helper (optional)
function toBinImage(x) { return ee.Image(x).gt(0); }

function exportMask(img, name){
var im = toBinImage(img).unmask(0).toByte().rename('mask');
Export.image.toDrive({
image: im,
description: name + '_GEOTIFF',
fileNamePrefix: name,
region: WE_forestsGeom,
scale: 30,
maxPixels: 1e13
});
}

/***** EXPORT (FOREST MASK extent) *****/

var FOREST_REGION = ee.Geometry(WE_forestsGeom).bounds(1);

// A) Do you want a GeoTIFF with 0 outside the forest (easy for GIS overlay)?
function toMask01_zerosOutsideForest(img){
  return ee.Image(img)
    .gt(0)                       // 1 where true
    .updateMask(WE_forestMask30) // keep forest only
    .unmask(0)                   // 0 outside forest (within the exported bbox)
    .toByte()
    .rename('mask');
}

// B) (alternative) Do you want NoData outside the forest (not 0)?
function toMask01_nodataOutsideForest(img){
  return ee.Image(img)
    .gt(0)
    .updateMask(WE_forestMask30)
    .toByte()
    .rename('mask');
}


var EXPORT_MODE = 'ZEROS'; // 'ZEROS' ή 'NODATA'

function exportForestMask(img, name){
  var out = (EXPORT_MODE === 'NODATA')
    ? toMask01_nodataOutsideForest(img)
    : toMask01_zerosOutsideForest(img);

  Export.image.toDrive({
    image: out,
    description: name + '_FOREST_GEOTIFF',
    fileNamePrefix: name,
    region: FOREST_REGION,   // ✅ forest extent (
    scale: 30,
    maxPixels: 1e13,
    fileFormat: 'GeoTIFF'
  });
}

// 6 exports
exportForestMask(A1_like_STRICT, 'A1_like_STRICT');
exportForestMask(A1_like_RELAX,  'A1_like_RELAX');

exportForestMask(A2_like_STRICT, 'A2_like_STRICT');
exportForestMask(A2_like_RELAX,  'A2_like_RELAX');

exportForestMask(A3_like_STRICT, 'A3_like_STRICT');
exportForestMask(A3_like_RELAX,  'A3_like_RELAX');


/***** PRECISION — VALIDATION — RECALL (+ MAJORITY DIAGNOSTICS)  *****/

/*** ======= COMMON HELPERS ======= ***/
function _bin(x){ return ee.Image(x).gt(0).unmask(0).toByte(); }
function _toBinBand(x){ return _bin(x).rename('b'); }

function _countOnes(maskImg, geom){
  return _bin(maskImg).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: geom,
    scale: 90,
    bestEffort: true,
    tileScale: 4,
    maxPixels: 2e9
  }).values().get(0);
}

function _countForest(geom){
  return ee.Image(1).updateMask(WE_forestMask30).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geom,
    scale: 90,
    bestEffort: true,
    tileScale: 4,
    maxPixels: 2e9
  }).values().get(0);
}

// % coverage inside the ROI (over the forested area)
function coveragePctROI(maskImg, roi){
  var geom = ee.Geometry(roi).intersection(WE_forestsGeom, 1);
  var sum = _countOnes(maskImg.updateMask(WE_forestMask30), geom);
  var tot = _countForest(geom);
  return ee.Number(sum).divide(ee.Number(tot)).multiply(100);
}

// Precision% = TP / predicted positives
// (TP always inside the ROI; predicted positives either inside the ROI or global)
var PRECISION_IN_ROI = true;


function precisionPct(predMask, truthMask, roi){
  var roiGeom = ee.Geometry(roi).intersection(WE_forestsGeom, 1);
  var scopeForPred = ee.Geometry(ee.Algorithms.If(PRECISION_IN_ROI, roiGeom, WE_forestsGeom));

  var predInScope = _countOnes(predMask.updateMask(WE_forestMask30), scopeForPred);
  var TP = _countOnes(
    _bin(predMask).and(_bin(truthMask)).updateMask(WE_forestMask30),
    roiGeom
  );

  return ee.Algorithms.If(
    ee.Number(predInScope).gt(0),
    ee.Number(TP).divide(ee.Number(predInScope)).multiply(100),
    null
  );
}

// Validation: % of the ROI area (over forest) that is covered by the mask
function pctAreaInROI(maskImg, roi){
  var pixArea = ee.Image.pixelArea().rename('area');
  var roiGeom = ee.Geometry(roi).intersection(WE_forestsGeom, 1);

  var hitArea = pixArea
    .updateMask(_toBinBand(maskImg))
    .updateMask(WE_forestMask30)
    .updateMask(ee.Image().paint(roiGeom, 1))
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: roiGeom,
      scale: 30,
      bestEffort: true,
      tileScale: 4,
      maxPixels: 2e9
    })
    .get('area');

  var denArea = pixArea
    .updateMask(WE_forestMask30)
    .updateMask(ee.Image().paint(roiGeom, 1))
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: roiGeom,
      scale: 30,
      bestEffort: true,
      tileScale: 4,
      maxPixels: 2e9
    })
    .get('area');

  return ee.Number(hitArea).divide(ee.Number(denArea)).multiply(100);
}

// Recall% = TP / total truth (within the ROI)
function recallPct(predMask, truthMask, roi){
  var roiGeom = ee.Geometry(roi).intersection(WE_forestsGeom, 1);

  var TP = _countOnes(
    _bin(predMask).and(_bin(truthMask)).updateMask(WE_forestMask30),
    roiGeom
  );

  var totalTruth = _countOnes(
    _bin(truthMask).updateMask(WE_forestMask30),
    roiGeom
  );

  return ee.Algorithms.If(
    ee.Number(totalTruth).gt(0),
    ee.Number(TP).divide(ee.Number(totalTruth)).multiply(100),
    null
  );
}


/*** ======= METRICS (STRICT/RELAX) per area ======= ***/
function metricsFeature(areaLabel, mode, likeMask, mortMask, roi){
  var cov = coveragePctROI(likeMask, roi);
  var pre = precisionPct(likeMask, mortMask, roi);
  var rec = recallPct(likeMask, mortMask, roi);
  var val = pctAreaInROI(likeMask, roi);

  return ee.Feature(null, {
    area: areaLabel,
    mode: mode,
    precision_scope: (PRECISION_IN_ROI ? 'ROI' : 'GLOBAL_FORESTS'),
    coverage_pct: cov,
    precision_pct: pre,
    recall_pct: rec,
    validation_pct: val
  });
}

var METRICS = ee.FeatureCollection([
  metricsFeature('Area1', 'STRICT', A1_like_STRICT, A1_mortMask, gArea1),
  metricsFeature('Area1', 'RELAX',  A1_like_RELAX,  A1_mortMask, gArea1),

  metricsFeature('Area2', 'STRICT', A2_like_STRICT, A2_mortMask, gArea2),
  metricsFeature('Area2', 'RELAX',  A2_like_RELAX,  A2_mortMask, gArea2),

  metricsFeature('Area3', 'STRICT', A3_like_STRICT, A3_mortMask, gArea3),
  metricsFeature('Area3', 'RELAX',  A3_like_RELAX,  A3_mortMask, gArea3)
]);

// STRICT & ένα RELAX print per area
print('Area1 — STRICT metrics', METRICS.filter(ee.Filter.and(
  ee.Filter.eq('area','Area1'),
  ee.Filter.eq('mode','STRICT')
)));
print('Area1 — RELAX metrics', METRICS.filter(ee.Filter.and(
  ee.Filter.eq('area','Area1'),
  ee.Filter.eq('mode','RELAX')
)));

print('Area2 — STRICT metrics', METRICS.filter(ee.Filter.and(
  ee.Filter.eq('area','Area2'),
  ee.Filter.eq('mode','STRICT')
)));
print('Area2 — RELAX metrics', METRICS.filter(ee.Filter.and(
  ee.Filter.eq('area','Area2'),
  ee.Filter.eq('mode','RELAX')
)));

print('Area3 — STRICT metrics', METRICS.filter(ee.Filter.and(
  ee.Filter.eq('area','Area3'),
  ee.Filter.eq('mode','STRICT')
)));
print('Area3 — RELAX metrics', METRICS.filter(ee.Filter.and(
  ee.Filter.eq('area','Area3'),
  ee.Filter.eq('mode','RELAX')
)));


/*** ======= COVERAGE & PRECISION SUMMARY TABLE ======= ***/
function reportPR(areaLabel, roi, likeStrict, likeRelax, mortMask){
  var cov_s = coveragePctROI(likeStrict, roi);
  var cov_r = coveragePctROI(likeRelax,  roi);

  var pre_s = precisionPct(likeStrict, mortMask, roi);
  var pre_r = precisionPct(likeRelax,  mortMask, roi);

  return ee.Feature(null, {
    area: areaLabel,
    precision_scope: (PRECISION_IN_ROI ? 'ROI' : 'GLOBAL_FORESTS'),
    coverage_STRICT_pct:  cov_s,
    precision_STRICT_pct: pre_s,
    coverage_RELAX_pct:   cov_r,
    precision_RELAX_pct:  pre_r
  });
}

var PR_TABLE = ee.FeatureCollection([
  reportPR('Area1', gArea1, A1_like_STRICT, A1_like_RELAX, A1_mortMask),
  reportPR('Area2', gArea2, A2_like_STRICT, A2_like_RELAX, A2_mortMask),
  reportPR('Area3', gArea3, A3_like_STRICT, A3_like_RELAX, A3_mortMask)
]);

print('Coverage & Precision — per area (STRICT/RELAX)', PR_TABLE);

var PR_TABLE_ROUND = PR_TABLE.map(function(f){
  var pS = f.get('precision_STRICT_pct');
  var pR = f.get('precision_RELAX_pct');

  var pS_round = ee.Algorithms.If(ee.Algorithms.IsEqual(pS, null), null, ee.Number(pS).format('%.2f'));
  var pR_round = ee.Algorithms.If(ee.Algorithms.IsEqual(pR, null), null, ee.Number(pR).format('%.2f'));

  return f.set({
    precision_STRICT_round: pS_round,
    precision_RELAX_round:  pR_round
  });
});

var precisionTable = ui.Chart.feature.byFeature(
    PR_TABLE_ROUND,
    'area',
    ['precision_STRICT_round', 'precision_RELAX_round']
  )
  .setChartType('Table')
  .setOptions({
    title: 'PRECISION — STRICT / RELAX (forest-only)',
    allowHtml: true
  });

print(precisionTable);


/*** ======= VALIDATION TABLE ======= ***/
var VALIDATION = ee.FeatureCollection([
  ee.Feature(null, {
    area: 'Area1',
    strict_pct: pctAreaInROI(A1_like_STRICT, gArea1),
    relax_pct:  pctAreaInROI(A1_like_RELAX,  gArea1)
  }),
  ee.Feature(null, {
    area: 'Area2',
    strict_pct: pctAreaInROI(A2_like_STRICT, gArea2),
    relax_pct:  pctAreaInROI(A2_like_RELAX,  gArea2)
  }),
  ee.Feature(null, {
    area: 'Area3',
    strict_pct: pctAreaInROI(A3_like_STRICT, gArea3),
    relax_pct:  pctAreaInROI(A3_like_RELAX,  gArea3)
  })
]);

print('VALIDATION — % κάλυψης μέσα στο AOI (forest-only)', VALIDATION);

var VALIDATION_ROUND = VALIDATION.map(function(f){
  return f.set({
    strict_pct_round: ee.Number(f.get('strict_pct')).format('%.2f'),
    relax_pct_round:  ee.Number(f.get('relax_pct')).format('%.2f')
  });
});

var validationTable = ui.Chart.feature.byFeature(
    VALIDATION_ROUND,
    'area',
    ['strict_pct_round', 'relax_pct_round']
  )
  .setChartType('Table')
  .setOptions({
    title: 'VALIDATION — % κάλυψης μέσα στο AOI (forest-only)',
    allowHtml: true
  });

print(validationTable);


/*** ======= RECALL TABLE ======= ***/
var RECALL = ee.FeatureCollection([
  ee.Feature(null, {
    area: 'Area1',
    recall_STRICT_pct: recallPct(A1_like_STRICT, A1_mortMask, gArea1),
    recall_RELAX_pct:  recallPct(A1_like_RELAX,  A1_mortMask, gArea1)
  }),
  ee.Feature(null, {
    area: 'Area2',
    recall_STRICT_pct: recallPct(A2_like_STRICT, A2_mortMask, gArea2),
    recall_RELAX_pct:  recallPct(A2_like_RELAX,  A2_mortMask, gArea2)
  }),
  ee.Feature(null, {
    area: 'Area3',
    recall_STRICT_pct: recallPct(A3_like_STRICT, A3_mortMask, gArea3),
    recall_RELAX_pct:  recallPct(A3_like_RELAX,  A3_mortMask, gArea3)
  })
]);

print('RECALL — % captured mortality (forest-only, within the AOI)', RECALL);

var RECALL_ROUND = RECALL.map(function(f){
  var rS = f.get('recall_STRICT_pct');
  var rR = f.get('recall_RELAX_pct');

  var rS_round = ee.Algorithms.If(ee.Algorithms.IsEqual(rS, null), null, ee.Number(rS).format('%.2f'));
  var rR_round = ee.Algorithms.If(ee.Algorithms.IsEqual(rR, null), null, ee.Number(rR).format('%.2f'));

  return f.set({
    recall_STRICT_round: rS_round,
    recall_RELAX_round:  rR_round
  });
});

var recallTable = ui.Chart.feature.byFeature(
    RECALL_ROUND,
    'area',
    ['recall_STRICT_round', 'recall_RELAX_round']
  )
  .setChartType('Table')
  .setOptions({
    title: 'RECALL — % of the confirmed mortality that is covered (STRICT / RELAX)',
    allowHtml: true
  });

print(recallTable);


/*** ======= MAJORITY DIAGNOSTICS (RELAX) ======= ***/
function getMortMaskByArea(areaLabel){
  if (areaLabel === 'Area1') return A1_mortMask;
  if (areaLabel === 'Area2') return A2_mortMask;
  return A3_mortMask;
}

function majorityDiagForArea(areaLabel, roi, likeRelaxMask){
  var roiGeom = ee.Geometry(roi).intersection(WE_forestsGeom, 1);

  // want helpers/vars from RELAX section:
  // RELAX_K, Q_RELAX_LOW/HIGH, RELAX_RULE
  // getQuantileRange, _clamp, _wrapDeg, _num, _bandRangeMask
  // WE_aspectWindowDeg, WE_slopeDeg, WE_TPI300, WE_TPI1000, WE_curv8
  // USE_TOPO, USE_ASPECT, USE_SLOPE, USE_TPI300, USE_TPI1000, USE_CURV8

  var topoD = ee.Dictionary(
    areaLabel === 'Area1' ? topoArea1.first().toDictionary()
    : areaLabel === 'Area2' ? topoArea2.first().toDictionary()
                            : topoArea3.first().toDictionary()
  );

  var K    = RELAX_K;
  var qLow = Q_RELAX_LOW;
  var qHigh= Q_RELAX_HIGH;

  var aspect_mu = _wrapDeg(_num(topoD,'aspect_circ_mean_deg',180));
  var aspect_sd = _num(topoD,'aspect_circ_std_deg',30).multiply(K);
  var m_aspect  = WE_aspectWindowDeg(aspect_mu, aspect_sd);

  var slopeQ = getQuantileRange(WE_slopeDeg, roiGeom, qLow, qHigh);
  var slope_lo = _clamp(slopeQ.lo, 0, 90);
  var slope_hi = _clamp(slopeQ.hi, 0, 90);
  var m_slope  = WE_slopeDeg.gte(slope_lo).and(WE_slopeDeg.lte(slope_hi));

  var tpi3Q = getQuantileRange(WE_TPI300, roiGeom, qLow, qHigh);
  var m_tpi3 = _bandRangeMask(WE_TPI300, tpi3Q.lo, tpi3Q.hi);

  var tpi1Q = getQuantileRange(WE_TPI1000, roiGeom, qLow, qHigh);
  var m_tpi1 = _bandRangeMask(WE_TPI1000, tpi1Q.lo, tpi1Q.hi);

  var curvQ = getQuantileRange(WE_curv8, roiGeom, qLow, qHigh);
  var m_curv8 = _bandRangeMask(WE_curv8, curvQ.lo, curvQ.hi);

  var mortMask = getMortMaskByArea(areaLabel);
  var likeCov  = coveragePctROI(likeRelaxMask, roi);

  function makeRow(name, useFlag, mFeat){
    if (!useFlag) return null;

    var compMask = mFeat.and(likeRelaxMask).updateMask(WE_forestMask30);

    var cov = coveragePctROI(compMask, roi);
    var share = ee.Algorithms.If(
      ee.Number(likeCov).gt(0),
      ee.Number(cov).divide(likeCov).multiply(100),
      null
    );

    var prec = precisionPct(compMask, mortMask, roi);
    var rec  = recallPct(compMask,  mortMask, roi);
    var val  = pctAreaInROI(compMask, roi);

    return ee.Feature(null, {
      area: areaLabel,
      mode: 'RELAX',
      majority_rule: RELAX_RULE,
      feature: name,
      coverage_pct: cov,
      share_of_relax_pct: share,
      precision_pct: prec,
      recall_pct: rec,
      validation_pct: val
    });
  }

  var rows = [];

  if (USE_TOPO && USE_ASPECT)  rows.push(makeRow('ASPECT',  true, m_aspect));
  if (USE_TOPO && USE_SLOPE)   rows.push(makeRow('SLOPE',   true, m_slope));
  if (USE_TOPO && USE_TPI300)  rows.push(makeRow('TPI_300', true, m_tpi3));
  if (USE_TOPO && USE_TPI1000) rows.push(makeRow('TPI_1000',true, m_tpi1));
  if (USE_TOPO && USE_CURV8)   rows.push(makeRow('CURV8',   true, m_curv8));

  rows = rows.filter(function(r){ return r !== null; });

  return ee.FeatureCollection(rows);
}

// diagnostics for the 3 aoi (RELAX)
var MAJ_DIAG_A1 = majorityDiagForArea('Area1', gArea1, A1_like_RELAX);
var MAJ_DIAG_A2 = majorityDiagForArea('Area2', gArea2, A2_like_RELAX);
var MAJ_DIAG_A3 = majorityDiagForArea('Area3', gArea3, A3_like_RELAX);

print('MAJORITY — per-feature diagnostics (Area1 — RELAX)', MAJ_DIAG_A1);
print('MAJORITY — per-feature diagnostics (Area2 — RELAX)', MAJ_DIAG_A2);
print('MAJORITY — per-feature diagnostics (Area3 — RELAX)', MAJ_DIAG_A3);

var majTable_A1 = ui.Chart.feature.byFeature(
    MAJ_DIAG_A1,
    'feature',
    ['coverage_pct','share_of_relax_pct','precision_pct','recall_pct','validation_pct']
  )
  .setChartType('Table')
  .setOptions({ title: 'MAJORITY — Area1 (RELAX, per characteristic)', allowHtml: true });
print(majTable_A1);

var majTable_A2 = ui.Chart.feature.byFeature(
    MAJ_DIAG_A2,
    'feature',
    ['coverage_pct','share_of_relax_pct','precision_pct','recall_pct','validation_pct']
  )
  .setChartType('Table')
  .setOptions({ title: 'MAJORITY — Area2 (RELAX, per characteristic)', allowHtml: true });
print(majTable_A2);

var majTable_A3 = ui.Chart.feature.byFeature(
    MAJ_DIAG_A3,
    'feature',
    ['coverage_pct','share_of_relax_pct','precision_pct','recall_pct','validation_pct']
  )
  .setChartType('Table')
  .setOptions({ title: 'MAJORITY — Area3 (RELAX, per characteristic)', allowHtml: true });
print(majTable_A3);

/***** RELAX — SIMPLE PRINT (per area) *****/

function printRelaxParams(areaLabel){
  var roiGeom = (areaLabel === 'Area1') ? gArea1
             : (areaLabel === 'Area2') ? gArea2
                                       : gArea3;

  var topoD = ee.Dictionary(
    areaLabel==='Area1' ? topoArea1.first().toDictionary()
  : areaLabel==='Area2' ? topoArea2.first().toDictionary()
                        : topoArea3.first().toDictionary()
  );

  var solD  = ee.Dictionary(
    areaLabel==='Area1' ? solArea1.first().toDictionary()
  : areaLabel==='Area2' ? solArea2.first().toDictionary()
                        : solArea3.first().toDictionary()
  );

  // EVENT thresholds (RELAX) από AUTO truth
  var thr_drop_nbr_relax  = ee.Number(WE_MIN_DROP_NBR).multiply(RELAX_DROP_FACTOR);
  var thr_drop_ndvi_relax = ee.Number(WE_MIN_DROP_NDVI).multiply(RELAX_DROP_FACTOR);

  // PERSIST thresholds (RELAX) όπως στο makeLikeMask (final)
  var thr_p1_truth_nbr_relax  = ee.Number(WE_P1_DELTA_NBR).multiply(RELAX_PERSIST_FACTOR);
  var thr_p1_truth_ndvi_relax = ee.Number(WE_P1_DELTA_NDVI).multiply(RELAX_PERSIST_FACTOR);

  var thr_frac_nbr  = ee.Number(WE_EVT_DELTA_NBR).multiply(PERSIST_DELTA_FRAC);
  var thr_frac_ndvi = ee.Number(WE_EVT_DELTA_NDVI).multiply(PERSIST_DELTA_FRAC);

  var thr_p1_final_nbr  = ee.Number(thr_p1_truth_nbr_relax).max(thr_frac_nbr).max(0);
  var thr_p1_final_ndvi = ee.Number(thr_p1_truth_ndvi_relax).max(thr_frac_ndvi).max(0);

  // RELAX quantile ranges (TOPO)
  var qLow  = Q_RELAX_LOW;
  var qHigh = Q_RELAX_HIGH;
  var K     = RELAX_K;

  var aspect_mu = _wrapDeg(_num(topoD,'aspect_circ_mean_deg',180));
  var aspect_hw = _num(topoD,'aspect_circ_std_deg',30).multiply(K);

  var slopeQ  = getQuantileRange(WE_slopeDeg, roiGeom, qLow, qHigh);
  var slopeLo = _clamp(slopeQ.lo, 0, 90);
  var slopeHi = _clamp(slopeQ.hi, 0, 90);

  var tpi3Q = getQuantileRange(WE_TPI300,  roiGeom, qLow, qHigh);
  var tpi1Q = getQuantileRange(WE_TPI1000, roiGeom, qLow, qHigh);
  var curvQ = getQuantileRange(WE_curv8,   roiGeom, qLow, qHigh);

  // SOLAR (όπως στο makeLikeMask)
  var solar_mu  = _num(solD,'mean',0.90);
  var solar_sd  = _num(solD,'stdDev',0.03);
  var solar_min = ee.Number(solar_mu).subtract(solar_sd).max(0).min(1);
  var useSolarThisArea = (areaLabel!=='Area2') && USE_SOLAR;
  var solar_min_used   = ee.Algorithms.If(useSolarThisArea, solar_min, null);

  // SPEI settings 
  var spei_mode = SPEI_MODE;
  var spei_thr  = WE_SPEI_THR_ALL;

  // GEOLOGY top3 codes
  var geocodes = geologyTopCodes(areaLabel, 3);

  // ELEV 
  var elevLo = null, elevHi = null;
  if (USE_ELEVATION){
    var elevObj = elevationMaskFromStats(areaLabel);
    elevLo = elevObj.lo;
    elevHi = elevObj.hi;
  }

  var out = ee.Dictionary({
    area: areaLabel,
    mode: 'RELAX',

    relax_rule: RELAX_RULE,
    relax_K: RELAX_K,
    q_low: Q_RELAX_LOW,
    q_high: Q_RELAX_HIGH,
    relax_drop_factor: RELAX_DROP_FACTOR,
    relax_persist_factor: RELAX_PERSIST_FACTOR,
    require_persist_post1: REQUIRE_PERSIST_POST1,
    persist_delta_frac: PERSIST_DELTA_FRAC,

    // thresholds used
    thr_drop_nbr_pct: thr_drop_nbr_relax,
    thr_drop_ndvi_pct: thr_drop_ndvi_relax,
    thr_p1_nbr_delta_final: thr_p1_final_nbr,
    thr_p1_ndvi_delta_final: thr_p1_final_ndvi,

    // topo ranges used
    aspect_mu_deg: aspect_mu,
    aspect_halfWidth_deg: aspect_hw,
    slope_lo_deg: slopeLo,
    slope_hi_deg: slopeHi,
    tpi300_lo: tpi3Q.lo,
    tpi300_hi: tpi3Q.hi,
    tpi1000_lo: tpi1Q.lo,
    tpi1000_hi: tpi1Q.hi,
    curv8_lo: curvQ.lo,
    curv8_hi: curvQ.hi,

    // solar/spei/geol/elev
    solar_min_used: solar_min_used,
    spei_mode: spei_mode,
    spei_thr: spei_thr,
    geol_top3_codes: geocodes,
    elev_lo: elevLo,
    elev_hi: elevHi
  });

  print('RELAX PARAMS — ' + areaLabel, out);
}

// run
printRelaxParams('Area1');
printRelaxParams('Area2');
printRelaxParams('Area3');


/***** =========================================================
 * HEATMAP (LOCKED / ZOOM-INDEPENDENT) + HEAT RISK + COVERAGE/RECALL
 * Fixes:
 *  - low values become TRANSPARENT (no “black carpet”)
 *  - palette matches your screenshot (no black)
 * ========================================================= *****/

// -----------------------
// 0) HARD-BIN helpers (0/1 unmasked)
// -----------------------
function _binHard(img){
  return ee.Image(img).gt(0).unmask(0).toByte();
}

// -----------------------
// 1) LOCKED GRID (stable projection + scale)
// -----------------------
var HEAT_SCALE = 90;
var HEAT_PROJ  = ee.Projection('EPSG:32636').atScale(HEAT_SCALE);

function _forestMaskLocked(){
  return ee.Image(WE_forestMask30)
    .unmask(0)
    .toByte()
    .reproject({crs: HEAT_PROJ});
}

function _countSumLocked(bin01, geom){
  return ee.Image(bin01)
    .reproject({crs: HEAT_PROJ})
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: geom,
      scale: HEAT_SCALE,
      bestEffort: true,
      tileScale: 4,
      maxPixels: 2e9
    }).values().get(0);
}

function _countForestLocked(geom){
  var f = _forestMaskLocked();
  return ee.Image(1).updateMask(f).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: geom,
    scale: HEAT_SCALE,
    bestEffort: true,
    tileScale: 4,
    maxPixels: 2e9
  }).values().get(0);
}

// -----------------------
// 2) GLOBAL RELAX (OR) — unify band names
// -----------------------
var GLOBAL_RELAX_B = ee.ImageCollection([
  _binHard(A1_like_RELAX).rename('b'),
  _binHard(A2_like_RELAX).rename('b'),
  _binHard(A3_like_RELAX).rename('b')
]).max().rename('b');

// -----------------------
// 3) HEAT params
// -----------------------

// ✅ Palette like your screenshot (low risk -> high risk)
// low risk  -------------------------------> high risk
var HEAT_PALETTE = [
  '2b2b6e', // deep indigo / dark (low)
  '5a2ca0', // purple
  'ff007f', // magenta
  'ff2b2b', // red
  'ff7a00', // orange
  'ffb000', // amber
  'ffff00'  // yellow (high)
];


// IMPORTANT: this is now a REAL mask threshold (low values become transparent)
var HEAT_MIN_SHOW = 0.20;   // tweak: 0.01..0.08 depending on taste

var HEAT_GAMMA    = 0.70;   // boosts faint signal a bit (optional)
var HEAT_DILATE_M = 120;    // optional “blob” expansion before smoothing

// ArcGIS-like stretch (percentiles)
var HEAT_STRETCH_PLO = 2;
var HEAT_STRETCH_PHI = 98;

// Gaussian kernel (constant radius)
var MAX_KERNEL_PX = 512;
var BASE_SCALE_M  = 30;

var RADIUS_M = 2000;
var SIGMA_M  = RADIUS_M / 3;

// safety cap
var maxRadiusM = MAX_KERNEL_PX * BASE_SCALE_M;
RADIUS_M = Math.min(RADIUS_M, maxRadiusM);
SIGMA_M  = Math.min(SIGMA_M,  maxRadiusM);

// -----------------------
// 4) HEATMAP — LOCKED
// -----------------------
var forested = GLOBAL_RELAX_B
  .updateMask(WE_forestMask30)
  .unmask(0)
  .toFloat()
  .reproject({crs: HEAT_PROJ});

if (HEAT_DILATE_M && HEAT_DILATE_M > 0){
  forested = forested
    .focal_max({radius: HEAT_DILATE_M, units: 'meters'})
    .reproject({crs: HEAT_PROJ})
    .updateMask(WE_forestMask30);
}

var k = ee.Kernel.gaussian({
  radius: RADIUS_M,
  sigma:  SIGMA_M,
  units: 'meters',
  normalize: true
});

var smooth = forested.convolve(k).reproject({crs: HEAT_PROJ});

// percentile stretch
var p = smooth.reduceRegion({
  reducer: ee.Reducer.percentile([HEAT_STRETCH_PLO, HEAT_STRETCH_PHI]),
  geometry: WE_forestsGeom,
  scale: 300,
  bestEffort: true,
  tileScale: 4,
  maxPixels: 2e9
});

var pLo = ee.Number(p.get('b_p' + HEAT_STRETCH_PLO));
var pHi = ee.Number(p.get('b_p' + HEAT_STRETCH_PHI)).max(pLo.add(1e-6));

// 0..1 stretched + gamma
var HEAT = smooth
  .unitScale(pLo, pHi)
  .clamp(0, 1)
  .pow(HEAT_GAMMA)
  .updateMask(WE_forestMask30)
  .rename('HEAT');

// ✅ KEY FIX: low values become TRANSPARENT (no “black” low-risk carpet)
HEAT = HEAT.updateMask(HEAT.gte(HEAT_MIN_SHOW));

Map.addLayer(
  HEAT,
  {min: 0, max: 1, palette: HEAT_PALETTE},
  'HEATMAP (LOCKED + transparent low)',
  true
);

// -----------------------
// 5) PRED mask from HEAT (for metrics)
// -----------------------
var HEAT_PRED = HEAT.gt(0).rename('pred');

// -----------------------
// 6) COVERAGE / RECALL (LOCKED, forest-only)
// -----------------------
function coveragePctLocked(predMask01, roiGeom){
  var g = ee.Geometry(roiGeom).intersection(WE_forestsGeom, 1);
  var f = _forestMaskLocked();

  var pred = _binHard(predMask01).reproject({crs: HEAT_PROJ}).updateMask(f);

  var ones = _countSumLocked(pred, g);
  var tot  = _countForestLocked(g);

  return ee.Number(ones).divide(ee.Number(tot)).multiply(100).clamp(0, 100);
}

function recallPctLocked(predMask01, truthMask01, roiGeom){
  var g = ee.Geometry(roiGeom).intersection(WE_forestsGeom, 1);
  var f = _forestMaskLocked();

  var pred  = _binHard(predMask01).reproject({crs: HEAT_PROJ}).updateMask(f);
  var truth = _binHard(truthMask01).reproject({crs: HEAT_PROJ}).updateMask(f);

  var TP = _countSumLocked(pred.and(truth), g);
  var TT = _countSumLocked(truth, g);

  return ee.Algorithms.If(
    ee.Number(TT).gt(0),
    ee.Number(TP).divide(ee.Number(TT)).multiply(100).clamp(0, 100),
    null
  );
}

// -----------------------
// 7) HEAT RISK (LOCKED)
// -----------------------
function riskPct_meanHeatLocked(polyGeom){
  var g = ee.Geometry(polyGeom).intersection(WE_forestsGeom, 1);
  var f = _forestMaskLocked();

  var heatLock = HEAT
    .unmask(0)
    .reproject({crs: HEAT_PROJ})
    .updateMask(f);

  var meanHeat = heatLock.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: g,
    scale: HEAT_SCALE,
    bestEffort: true,
    tileScale: 4,
    maxPixels: 2e9
  }).get('HEAT');

  return ee.Number(meanHeat).multiply(100).clamp(0, 100);
}

function riskPct_hotShareLocked(polyGeom, thr){
  var g = ee.Geometry(polyGeom).intersection(WE_forestsGeom, 1);
  var f = _forestMaskLocked();

  var hot = HEAT
    .unmask(0)
    .gte(thr)
    .toByte()
    .reproject({crs: HEAT_PROJ})
    .updateMask(f);

  var ones = _countSumLocked(hot, g);
  var tot  = _countForestLocked(g);

  return ee.Number(ones).divide(ee.Number(tot)).multiply(100).clamp(0, 100);
}

// thresholds
var HEAT_T_RED  = 0.60;
var HEAT_T_PURP = 0.80;

// OPTIONAL: polyA / polyB metrics (if you have them)
function _truthMaskForPoly(polyGeom){
  return ee.Image(1).clip(polyGeom).updateMask(WE_forestMask30).rename('truth');
}

if (typeof polyA !== 'undefined'){
  var polyA_truth = _truthMaskForPoly(polyA);

  print('HEATMAP metrics — polyA (LOCKED)', ee.Dictionary({
    coverage_pct: coveragePctLocked(HEAT_PRED, polyA),
    recall_pct:   recallPctLocked(HEAT_PRED, polyA_truth, polyA)
  }));

  print('HEAT RISK — polyA (LOCKED)', ee.Dictionary({
    meanRiskPct:            riskPct_meanHeatLocked(polyA),
    hotspot_ge_0_60_pct:    riskPct_hotShareLocked(polyA, HEAT_T_RED),
    hotspot_ge_0_80_pct:    riskPct_hotShareLocked(polyA, HEAT_T_PURP)
  }));
}

if (typeof polyB !== 'undefined'){
  var polyB_truth = _truthMaskForPoly(polyB);

  print('HEATMAP metrics — polyB (LOCKED)', ee.Dictionary({
    coverage_pct: coveragePctLocked(HEAT_PRED, polyB),
    recall_pct:   recallPctLocked(HEAT_PRED, polyB_truth, polyB)
  }));

  print('HEAT RISK — polyB (LOCKED)', ee.Dictionary({
    meanRiskPct:            riskPct_meanHeatLocked(polyB),
    hotspot_ge_0_60_pct:    riskPct_hotShareLocked(polyB, HEAT_T_RED),
    hotspot_ge_0_80_pct:    riskPct_hotShareLocked(polyB, HEAT_T_PURP)
  }));
}


//Covered forest from heatmap

/***** % COVERAGE of ALL CyprusForestsAll by HEATMAP (forest-only) *****/

// Use the SAME locked projection/scale you already defined:
var HEAT_SCALE = 90;
var HEAT_PROJ  = ee.Projection('EPSG:32636').atScale(HEAT_SCALE);

function forestMaskLocked(){
  return ee.Image(WE_forestMask30)
    .unmask(0).toByte()
    .reproject({crs: HEAT_PROJ});
}

// Count forest pixels (denominator)
function countForestAll(){
  var f = forestMaskLocked();
  return ee.Image(1).updateMask(f).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: WE_forestsGeom,
    scale: HEAT_SCALE,
    bestEffort: true,
    tileScale: 4,
    maxPixels: 2e9
  }).values().get(0);
}

// Count covered forest pixels by a binary mask (numerator)
function countCoveredAll(binMask01){
  var f = forestMaskLocked();
  var b = ee.Image(binMask01).gt(0).unmask(0).toByte()
    .reproject({crs: HEAT_PROJ})
    .updateMask(f);

  return b.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: WE_forestsGeom,
    scale: HEAT_SCALE,
    bestEffort: true,
    tileScale: 4,
    maxPixels: 2e9
  }).values().get(0);
}

// ✅ Heat “presence” mask = where HEAT is visible (non-transparent)
var HEAT_VISIBLE = HEAT.mask().rename('heat_vis');   // 1 where HEAT exists, 0 elsewhere

var forestTot   = ee.Number(countForestAll());
var heatCovered = ee.Number(countCoveredAll(HEAT_VISIBLE));

var heatCoveragePct_allForest = heatCovered.divide(forestTot).multiply(100);

print('HEATMAP coverage of ALL CyprusForestsAll (forest-only) % =',
      heatCoveragePct_allForest.format('%.2f'));
      
/***** % of ALL CyprusForestsAll (forest-only) with HEAT >= thresholds *****/

var HEAT_SCALE = 90;
var HEAT_PROJ  = ee.Projection('EPSG:32636').atScale(HEAT_SCALE);

function forestMaskLocked(){
  return ee.Image(WE_forestMask30)
    .unmask(0).toByte()
    .reproject({crs: HEAT_PROJ});
}

function countForestAll(){
  var f = forestMaskLocked();
  return ee.Image(1).updateMask(f).reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: WE_forestsGeom,
    scale: HEAT_SCALE,
    bestEffort: true,
    tileScale: 4,
    maxPixels: 2e9
  }).values().get(0);
}

function pctForestHeatGE(thr){
  var f = forestMaskLocked();

  var hot = HEAT.unmask(0).gte(thr).toByte()
    .reproject({crs: HEAT_PROJ})
    .updateMask(f);

  var ones = ee.Number(hot.reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: WE_forestsGeom,
    scale: HEAT_SCALE,
    bestEffort: true,
    tileScale: 4,
    maxPixels: 2e9
  }).values().get(0));

  var tot = ee.Number(countForestAll());
  return ones.divide(tot).multiply(100);
}

print('% of ALL CyprusForestsAll (forest-only) with HEAT >= 0.60 =',
      pctForestHeatGE(0.60).format('%.2f'));

print('% of ALL CyprusForestsAll (forest-only) with HEAT >= 0.80 =',
      pctForestHeatGE(0.80).format('%.2f'));

