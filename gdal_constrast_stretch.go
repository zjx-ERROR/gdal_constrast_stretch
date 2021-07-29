package main

import (
	"log"
	"math"

	"reflect"

	"github.com/lukeroth/gdal"
)

type Binning struct {
	Nbins  int
	Offset float64
	Scale  float64
}

func (self *Binning) ToBin(v float64) int {
	if math.IsInf(v, -1) {
		return 0
	}
	if math.IsInf(v, 1) {
		return self.Nbins - 1
	}
	bin_dbl := math.Round((v - self.Offset) / self.Scale)
	if math.IsNaN(bin_dbl) {
		log.Fatal("nan in to_bin")
	}
	if bin_dbl < 0 {
		return 0
	}
	if bin_dbl > float64(self.Nbins)-1 {
		return self.Nbins - 1
	}
	return int(bin_dbl)
}

func (self *Binning) FromBin(i int) float64 {
	return float64(i)*self.Scale + self.Offset
}

type Histogram struct {
	Binning
	Max, Min, Mean, Stddev float64
	DataCount              uint
	NdvCount               uint
	Counts                 []uint
}

func ComputeMinmax(src_bands []gdal.RasterBand, ndv_def *NdvDef, w, h int) [][2]float64 {
	band_count := len(src_bands)
	minmax := make([][2]float64, band_count)
	blocksize_x_int, blocksize_y_int := src_bands[0].BlockSize()
	block_len := blocksize_x_int * blocksize_y_int
	buf_in := make([][]float64, band_count)
	for band_idx := 0; band_idx < band_count; band_idx++ {
		buf_in[band_idx] = make([]float64, block_len)
	}
	ndv_mask := make([]uint8, block_len)
	got_data := make([]bool, block_len)
	for boff_y := 0; boff_y < h; boff_y += blocksize_y_int {
		bsize_y := blocksize_y_int
		if bsize_y+boff_y > h {
			bsize_y = h - boff_y
		}
		for boff_x := 0; boff_x < w; boff_x += blocksize_x_int {
			bsize_x := blocksize_x_int
			if bsize_x+boff_x > w {
				bsize_x = w - boff_x
			}
			block_len = bsize_x * bsize_y

			for band_idx := 0; band_idx < band_count; band_idx++ {
				src_bands[band_idx].IO(gdal.Read, boff_x, boff_y, bsize_x, bsize_y, buf_in[band_idx], bsize_x, bsize_y, 0, 0)
			}
			ndv_def.GetNdvMaskC(buf_in, ndv_mask, block_len)
			for band_idx := 0; band_idx < band_count; band_idx++ {
				for i := 0; i < block_len; i++ {
					if ndv_mask[i] != 0 {
						continue
					}
					v := buf_in[band_idx][i]
					if math.IsNaN(v) || math.IsInf(v, 1) {
						continue
					}
					if !got_data[band_idx] {
						minmax[band_idx][0] = v
						minmax[band_idx][1] = v
						got_data[band_idx] = true
					}
					if v < minmax[band_idx][0] {
						minmax[band_idx][0] = v
					}
					if v > minmax[band_idx][1] {
						minmax[band_idx][1] = v
					}
				}
			}
		}
	}
	return minmax
}

func ComputeHistogram(src_bands []gdal.RasterBand, ndv_def *NdvDef, w, h int, binnings []Binning) []Histogram {
	band_count := len(src_bands)
	histograms := make([]Histogram, band_count)
	for band_idx := 0; band_idx < band_count; band_idx++ {
		histograms[band_idx].Binning = binnings[band_idx]
		histograms[band_idx].Counts = make([]uint, binnings[band_idx].Nbins)
	}

	blocksize_x_int, blocksize_y_int := src_bands[0].BlockSize()
	block_len := blocksize_x_int * blocksize_y_int
	buf_in := make([][]float64, band_count)
	for i := range buf_in {
		buf_in[i] = make([]float64, block_len)
	}
	ndv_mask := make([]uint8, block_len)
	first_valid_pixel := make([]bool, band_count)

	for boff_y := 0; boff_y < h; boff_y += blocksize_y_int {
		bsize_y := blocksize_y_int
		if bsize_y+boff_y > h {
			bsize_y = h - boff_y
		}
		for boff_x := 0; boff_x < w; boff_x += blocksize_x_int {
			bsize_x := blocksize_x_int
			if bsize_x+boff_x > w {
				bsize_x = w - boff_x
			}
			block_len = bsize_x * bsize_y

			for band_idx := 0; band_idx < band_count; band_idx++ {
				src_bands[band_idx].IO(gdal.Read, boff_x, boff_y, bsize_x, bsize_y, buf_in[band_idx], bsize_x, bsize_y, 0, 0)
			}
			ndv_def.GetNdvMaskC(buf_in, ndv_mask, block_len)

			for band_idx := 0; band_idx < band_count; band_idx++ {
				hg := &histograms[band_idx]
				p := buf_in[band_idx]
				for i := 0; i < block_len; i++ {
					if ndv_mask[i] != 0 {
						hg.NdvCount++
					} else {
						v := p[i]
						hg.Counts[hg.Binning.ToBin(v)]++
						if !first_valid_pixel[band_idx] {
							hg.Min = v
							hg.Max = v
							first_valid_pixel[band_idx] = true
						}
						if v < hg.Min {
							hg.Min = v
						}
						if v > hg.Max {
							hg.Max = v
						}
					}
				}
			}
		}

	}
	for band_idx := 0; band_idx < band_count; band_idx++ {
		hg := &histograms[band_idx]
		var accum float64
		for i := 0; i < hg.Binning.Nbins; i++ {
			cnt := hg.Counts[i]
			hg.DataCount += cnt
			accum += hg.Binning.FromBin(i) * float64(cnt)
		}
		hg.Mean = accum / float64(hg.DataCount)
		var var_accum float64
		for i := 0; i < hg.Binning.Nbins; i++ {
			cnt := hg.Counts[i]
			v := hg.Binning.FromBin(i)
			var_accum += (v - hg.Mean) * (v - hg.Mean) * float64(cnt)
		}
		hg.Stddev = math.Sqrt(var_accum / float64(hg.DataCount))
	}
	return histograms
}

func get_scale_from_percentile(histogram *Histogram, output_range int, from_percentile, to_percentile float64, scale_out, offset_out *float64) {
	start_count := uint(float64(histogram.DataCount) * from_percentile)
	end_count := uint(float64(histogram.DataCount) * to_percentile)
	var cnt uint = 0
	from_idx := -1
	to_idx := -1
	for i := 0; i < histogram.Binning.Nbins; i++ {
		if cnt <= start_count {
			from_idx = i
		}
		cnt += histogram.Counts[i]
		if cnt <= end_count {
			to_idx = i
		} else {
			break
		}
	}
	if from_idx < 0 || to_idx < 0 {
		log.Fatal("impossible: could not find window")
	}
	if from_idx == to_idx {
		from_idx = 0
		to_idx = histogram.Binning.Nbins - 1
	}

	from_val := histogram.Binning.FromBin(from_idx)
	to_val := histogram.Binning.FromBin(to_idx)

	*scale_out = float64(output_range-1)/to_val - from_val
	*offset_out = from_val
}

func invert_histogram(src_h_in *Histogram, dst_h []float64, output_range uint8) []uint8 {
	var (
		src_h       []uint
		pixel_count uint    = 0
		j           uint8   = 0
		src_total   float64 = 0
		dst_total   float64 = 0
	)
	copy(src_h_in.Counts, src_h)

	for i := 0; i < len(src_h); i++ {
		pixel_count += src_h[i]
	}
	out_h := make([]uint8, len(src_h))

	for i := 0; i < len(src_h); i++ {
		out_h[i] = j
		src_total += float64(src_h[i])
		for j < output_range-1 && dst_total < src_total {
			j++
			dst_total += dst_h[j] * float64(pixel_count)
		}
	}
	return out_h

}

func gen_gaussian(variance float64, bin_count int) []float64 {
	arr := make([]float64, bin_count)
	var total float64 = 0
	for i := 0; i < bin_count; i++ {
		if variance == 0.0 {
			arr[i] = 1
		} else {
			x := float64(i-bin_count/2) / variance
			arr[i] = math.Exp(-x * x)
		}
		total += arr[i]
	}
	for i := 0; i < bin_count; i++ {
		arr[i] /= total
	}
	return arr
}

func invert_histogram_to_gaussian(histogram_in *Histogram, variance float64, output_range int) []uint8 {
	gaussian := gen_gaussian(variance, output_range)
	return invert_histogram(histogram_in, gaussian, uint8(output_range))
}

func copyGeoCode(dst_ds, src_ds *gdal.Dataset) {
	affine := src_ds.GeoTransform()

	if len(affine) == 0 {
		affine = dst_ds.GeoTransform()
	}
	dst_ds.SetProjection(src_ds.Projection())
}

func usage() {

}

type Options struct {
	OutputFormat   string //-of
	Stddev         bool
	Percentile     bool
	Histeq         bool
	DumpHistogram  bool
	DstAvg         float64 //-linear-stretch
	DstStddev      float64
	FromPercentile float64 //-percentile-range
	ToPercentile   float64 //-percentile-range
	NdvLong        int64   //-outndv
	OutNdv         uint8   //-outndv
	SrcFn          string
	DstFn          string
	Ndv            [][2]float64
	ValidRange     [][2]float64
}

func (self *Options) handle() {
	var (
		ModeStddev        int = 0
		ModePercentile    int = 0
		ModeHisteq        int = 0
		ModeDumpHistogram int = 0
	)
	if self.Stddev {
		ModeStddev = 1
	} else {
		self.DstAvg = -1
		self.DstStddev = -1
	}
	if self.Percentile {
		ModePercentile = 1
	} else {
		self.FromPercentile = -1
		self.ToPercentile = -1
	}
	if self.Histeq {
		ModeHisteq = 1
	}
	if self.DumpHistogram {
		ModeDumpHistogram = 1
	}
	if ModeDumpHistogram+ModeHisteq+ModePercentile+ModeStddev > 1 || ModeDumpHistogram+ModeHisteq+ModePercentile+ModeStddev == 0 {
		log.Fatal("one mode to choose")
	}

	if len(self.SrcFn) == 0 {
		log.Fatal("missing srcfn")
	}
	if (len(self.DstFn) == 0) != self.DumpHistogram {
		log.Fatal("missing dstfn")
	}
	if self.Stddev && (self.DstAvg < 0 || self.DstStddev < 0) {
		log.Fatal("wrong agrs")
	}
	if self.Percentile && !(0 <= self.FromPercentile && self.FromPercentile < self.ToPercentile && self.ToPercentile <= 1) {
		log.Fatal("wrong args")
	}
	if len(self.OutputFormat) == 0 {
		self.OutputFormat = "GTiff"
	}
}

func Run(opt *Options) {
	opt.handle()
	ndv_def := NdvDef{}

	if len(opt.Ndv) > 0 && len(opt.ValidRange) > 0 {
		log.Fatal("you cannot use both Ndv and ValidRange")
	} else {
		ndv_def.Invert = len(opt.ValidRange) > 0
	}
	ndvslab := NdvSlab{}
	switch {
	case len(opt.Ndv) > 0:
		ndvslab.RangeByBand = opt.Ndv
	case len(opt.ValidRange) > 0:
		ndvslab.RangeByBand = opt.ValidRange
	}
	if !ndvslab.Empty() {
		ndv_def.Slabs = append(ndv_def.Slabs, ndvslab)
	}

	src_ds, err := gdal.Open(opt.SrcFn, gdal.ReadOnly)
	defer src_ds.Close()

	if err != nil {
		log.Fatal(err)
	}
	w := src_ds.RasterXSize()
	h := src_ds.RasterYSize()
	if w == 0 || h == 0 {
		log.Fatal("missing width/height")
	}
	src_band_count := src_ds.RasterCount()
	log.Printf("Input size is %d, %d, %d\n", w, h, src_band_count)

	var bandlist []int
	for i := 0; i < src_band_count; i++ {
		bandlist = append(bandlist, i+1)
	}
	dst_band_count := len(bandlist)

	if ndv_def.Empty() {
		var tmp [][2]float64
		band_count := src_ds.RasterCount()
		for _, v := range bandlist {
			if v < 1 || v > band_count {
				log.Fatal("bandid out of range")
			}
			band := src_ds.RasterBand(v)
			val, ok := band.NoDataValue()
			if ok {
				tmp = append(tmp, [2]float64{val, val})
			}
		}
		if len(tmp) > 0 {
			ndvslab.RangeByBand = tmp
			ndv_def.Slabs = append(ndv_def.Slabs, ndvslab)
		}
	}

	dst_driver, err := gdal.GetDriverByName(opt.OutputFormat)
	if err != nil {
		log.Fatal(err)
	}
	var src_bands, dst_bands []gdal.RasterBand
	for _, v := range bandlist {
		src_bands = append(src_bands, src_ds.RasterBand(v))
	}
	var binnings = make([]Binning, dst_band_count)
	var minmax [][2]float64
	for band_idx := 0; band_idx < dst_band_count; band_idx++ {
		var binning = &binnings[band_idx]
		dt := src_bands[band_idx].RasterDataType()
		switch dt {
		case gdal.Byte:
			binning.Nbins = 256
			binning.Offset = 0
			binning.Scale = 1
		case gdal.UInt16:
			binning.Nbins = 65536
			binning.Offset = 0
			binning.Scale = 1
		case gdal.Int16:
			binning.Nbins = 65536
			binning.Offset = -32768
			binning.Scale = 1
		default:
			if len(minmax) == 0 {
				minmax = ComputeMinmax(src_bands, &ndv_def, w, h)
			}
			binning.Nbins = 10000000
			binning.Offset = minmax[band_idx][0]
			binning.Scale = (minmax[band_idx][1] - minmax[band_idx][0]) / float64(binning.Nbins-1)
		}
	}
	print("\nComputing histogram...\n")

	histograms := ComputeHistogram(src_bands, &ndv_def, w, h, binnings)
	for band_idx := 0; band_idx < dst_band_count; band_idx++ {
		hg := &histograms[band_idx]
		log.Printf("band %d: min=%f, max=%f, mean=%f, stddev=%f, valid_count=%d, ndv_count=%d\n", band_idx+1, hg.Min, hg.Max, hg.Mean, hg.Stddev, hg.DataCount, hg.NdvCount)
		if opt.DumpHistogram {
			for i := 0; i < hg.Binning.Nbins; i++ {
				log.Printf("bin %d: val=%f cnt=%d\n", i, hg.Binning.FromBin(i), hg.Counts[i])
			}
		}
	}

	if opt.DumpHistogram {
		return
	}

	dst_ds := dst_driver.Create(opt.DstFn, w, h, dst_band_count, gdal.Byte, nil)
	defer dst_ds.Close()
	if reflect.DeepEqual(dst_ds, gdal.Dataset{}) {
		log.Fatal("couldn't create output")
	}
	copyGeoCode(&dst_ds, &src_ds)
	for band_idx := 0; band_idx < dst_band_count; band_idx++ {
		dst_bands = append(dst_bands, dst_ds.RasterBand(band_idx+1))
	}

	var (
		use_table    bool
		output_range = 256
		xform_table  = make([][]uint8, dst_band_count)
		lin_scales   = make([]float64, dst_band_count)
		lin_offsets  = make([]float64, dst_band_count)
	)
	if opt.Histeq {
		use_table = true
		for band_idx := 0; band_idx < dst_band_count; band_idx++ {
			xform_table[band_idx] = invert_histogram_to_gaussian(&histograms[band_idx], opt.DstStddev, output_range)
		}
	} else {
		use_table = false
		if opt.Percentile {
			for band_idx := 0; band_idx < dst_band_count; band_idx++ {
				get_scale_from_percentile(&histograms[band_idx], output_range, opt.FromPercentile, opt.ToPercentile, &lin_scales[band_idx], &lin_offsets[band_idx])
			}
		} else if opt.Stddev {
			for band_idx := 0; band_idx < dst_band_count; band_idx++ {
				hg := &histograms[band_idx]
				if hg.Stddev == 0 {
					lin_scales[band_idx] = 0
				} else {
					lin_scales[band_idx] = opt.DstStddev / hg.Stddev
				}
				lin_offsets[band_idx] = hg.Mean - opt.DstAvg/lin_scales[band_idx]
			}
		} else {
			print("\nWarning: no transformation was specified!  I'll just cast the input to 8-bit.\n")
			for band_idx := 0; band_idx < dst_band_count; band_idx++ {
				lin_scales[band_idx] = 1
				lin_offsets[band_idx] = 0
			}
		}

	}
	if use_table {
		if !ndv_def.Empty() {
			for j := range xform_table {
				for i := 0; i < len(xform_table[j]); i++ {
					va := xform_table[j][i]
					if va == opt.OutNdv {
						if opt.OutNdv < uint8(output_range)/2 {
							va++
						} else {
							va--
						}
					}
					xform_table[j][i] = va
				}
			}
		}
	}

	print("\nComputing output...\n")

	blocksize_x_int, blocksize_y_int := src_bands[0].BlockSize()
	block_len := blocksize_x_int * blocksize_y_int
	buf_in := make([][]float64, dst_band_count)
	buf_out := make([][]uint8, dst_band_count)
	for band_idx := 0; band_idx < dst_band_count; band_idx++ {
		buf_in[band_idx] = make([]float64, block_len)
		buf_out[band_idx] = make([]uint8, block_len)
	}
	ndv_mask := make([]uint8, block_len)

	for boff_y := 0; boff_y < h; boff_y += blocksize_y_int {
		bsize_y := blocksize_y_int
		if bsize_y+boff_y > h {
			bsize_y = h - boff_y
		}
		for boff_x := 0; boff_x < w; boff_x += blocksize_x_int {
			bsize_x := blocksize_x_int
			if bsize_x+boff_x > w {
				bsize_x = w - boff_x
			}
			block_len = bsize_x * bsize_y
			for band_idx := 0; band_idx < dst_band_count; band_idx++ {
				src_bands[band_idx].IO(gdal.Read, boff_x, boff_y, bsize_x, bsize_y, buf_in[band_idx], bsize_x, bsize_y, 0, 0)
			}

			ndv_def.GetNdvMaskC(buf_in, ndv_mask, block_len)

			for band_idx := 0; band_idx < dst_band_count; band_idx++ {
				p_in := &buf_in[band_idx]
				p_out := &buf_out[band_idx]
				p_ndv := &ndv_mask
				if use_table {
					xfrom := &xform_table[band_idx]
					binning := binnings[band_idx]
					for i := 0; i < block_len; i++ {
						if (*p_ndv)[i] != 0 {
							(*p_out)[i] = opt.OutNdv
						} else {
							(*p_out)[i] = (*xfrom)[binning.ToBin((*p_in)[i])]
						}
					}
				} else {
					scale := lin_scales[band_idx]
					offset := lin_offsets[band_idx]
					for i := 0; i < block_len; i++ {
						if (*p_ndv)[i] != 0 {
							(*p_out)[i] = opt.OutNdv
						} else {
							out_dbl := ((*p_in)[i] - offset) * scale
							var v uint8
							if out_dbl < 0 {
								v = 0
							} else if out_dbl > float64(output_range)-1 {
								v = uint8(output_range) - 1
							} else {
								v = uint8(out_dbl)
							}
							if v == opt.OutNdv {
								if opt.OutNdv < uint8(output_range)/2 {
									v++
								} else {
									v--
								}
							}
							(*p_out)[i] = v
						}
					}
				}
				dst_bands[band_idx].IO(gdal.Write, boff_x, boff_y, bsize_x, bsize_y, buf_out[band_idx], bsize_x, bsize_y, 0, 0)
			}
		}
	}

}
