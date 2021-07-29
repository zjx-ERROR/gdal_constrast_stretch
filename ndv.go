package main

import (
	"bytes"
	"encoding/binary"
	"fmt"
	"log"
	"unsafe"

	"github.com/lukeroth/gdal"
)

type NdvSlab struct {
	RangeByBand [][2]float64
}

func (self *NdvSlab) Empty() bool {
	return len(self.RangeByBand) == 0
}

func contains_templated(interval [2]float64, p interface{}) bool {
	switch p.(type) {
	case byte:
		return float64(p.(byte)) >= interval[0] && float64(p.(byte)) <= interval[1]
	case uint16:
		return float64(p.(uint16)) >= interval[0] && float64(p.(uint16)) <= interval[1]
	case int16:
		return float64(p.(int16)) >= interval[0] && float64(p.(int16)) <= interval[1]
	case int32:
		return float64(p.(int32)) >= interval[0] && float64(p.(int32)) <= interval[1]
	case uint32:
		return float64(p.(uint32)) >= interval[0] && float64(p.(uint32)) <= interval[1]
	case float32:
		return float64(p.(float32)) >= interval[0] && float64(p.(float32)) <= interval[1]
	case float64:
		return p.(float64) >= interval[0] && p.(float64) <= interval[1]
	case complex64:
		return float64(real(p.(complex64))) >= interval[0] && float64(real(p.(complex64))) <= interval[1]
	case complex128:
		return real(p.(complex128)) >= interval[0] && real(p.(complex128)) <= interval[1]
	default:
		log.Fatal("interval doesnot match any type")
		return false
	}
}

type NdvDef struct {
	Slabs  []NdvSlab
	Invert bool
}

func (self *NdvDef) PrintUsage() {
	print(
		"No-data values:\n",
		"  Ndv [[val val]]                                                          Set a no-data value\n",
		"  Ndv '[[val1 val1] [val2 val2] [val3 val3] ...'                          Set a no-data value using all input bands\n",
		"  Ndv '[[valMin1 valMax1] [valMin2 valMax2] [valMin3 valMax3] ...'         Set a range of no-data values\n",
		"  Ndv (-Inf and Inf are allowed; [[math.Inf(-1) math.Inf(1)]])\n",
		"  ValidRange '[[valMin1 valMax1] [valMin2 valMax2] [valMin3 valMax3] ...'  Set a range of valid data values\n",
	)
}

func (self *NdvDef) DebugPrint() {
	print("=== NDV\n")
	for i, iv := range self.Slabs {
		for j, jv := range iv.RangeByBand {
			fmt.Printf("range %d,%d = [%e,%e]\n", i, j, jv[0], jv[1])
		}
	}
	print("=== end NDV\n")
}

func (self *NdvDef) Empty() bool {
	return len(self.Slabs) == 0
}

func (self *NdvDef) IsInert() bool {
	return self.Invert
}

func (self *NdvDef) GetNdvMaskA(band float64, dt gdal.DataType, mask_out []uint8, num_pixels int) {
	var (
		bands   []float64
		dt_list []gdal.DataType
	)
	bands = append(bands, band)
	dt_list = append(dt_list, dt)
	self.GetNdvMaskB(bands, dt_list, mask_out, num_pixels)
}

func (self *NdvDef) GetNdvMaskB(bands []float64, dt_list []gdal.DataType, mask_out []uint8, num_pixels int) {
	var (
		in_p     [][]uint8
		dt_sizes []int
	)
	for i, v := range bands {
		in_p = append(in_p, convertToBytes(v))
		dt_sizes = append(dt_sizes, (dt_list[i] / 8).Size())
	}
	// MemSet(unsafe.Pointer(&mask_out), 0, uintptr(num_pixels))
	for pix_idx := 0; pix_idx < num_pixels; pix_idx++ {
		mask_out[pix_idx] = 0
		for j := range bands {
			if gdal_scalar_pointer_isnan(in_p[j][0], dt_list[j]) {
				mask_out[pix_idx] = 1
			}
		}
		for _, v := range self.Slabs {
			var all_match uint8 = 1
			for j := 0; j < len(bands); j++ {
				var k int
				if len(v.RangeByBand) == 1 {
					k = 0
				} else {
					k = j
				}
				if !contains_templated(v.RangeByBand[k], dt_list[j]) {
					all_match = 0
				}
			}
			mask_out[pix_idx] |= all_match
		}

		if self.Invert {
			if mask_out[pix_idx] == 0 {
				mask_out[pix_idx] = 1
			} else {
				mask_out[pix_idx] = 0
			}
		}

		for band_idx := 0; band_idx < len(bands); band_idx++ {
			in_p[band_idx] = in_p[band_idx][dt_sizes[band_idx]:]
		}
	}
}

func (self *NdvDef) GetNdvMaskC(bands [][]float64, mask_out []uint8, num_pixels int) {
	dt := gdal.Float64
	var (
		band_p  []float64
		dt_list []gdal.DataType
	)
	for _, v := range bands {
		band_p = append(band_p, v[0])
		dt_list = append(dt_list, dt)
	}
	self.GetNdvMaskB(band_p, dt_list, mask_out, num_pixels)
}

func (self *NdvDef) GetNdvMaskD(bands [][]float64, dt_list []gdal.DataType, mask_out []uint8, num_pixels int) {
	var band_p []float64
	for _, v := range bands {
		band_p = append(band_p, v[0])
	}
	self.GetNdvMaskB(band_p, dt_list, mask_out, num_pixels)
}

func gdal_scalar_pointer_isnan(p interface{}, dt gdal.DataType) bool {
	var ok bool
	switch dt {
	case gdal.Byte:
		_, ok = p.(uint8)
	case gdal.UInt16:
		_, ok = p.(uint16)
	case gdal.Int16:
		_, ok = p.(int16)
	case gdal.UInt32:
		_, ok = p.(uint32)
	case gdal.Int32:
		_, ok = p.(int32)
	case gdal.Float32:
		_, ok = p.(float32)
	case gdal.Float64:
		_, ok = p.(float64)
	// case gdal.CInt16:

	// case gdal.CInt32:
	case gdal.CFloat32:
		_, ok = p.(complex64)
	case gdal.CFloat64:
		_, ok = p.(complex128)
	default:
		log.Fatal("unrecognized datatype")
	}
	return ok
}

func convertToBytes(n interface{}) []uint8 {
	var buf bytes.Buffer
	err := binary.Write(&buf, binary.BigEndian, n)
	if err != nil {
		log.Panic(err)
	}
	return buf.Bytes()
}

func MemSet(s unsafe.Pointer, c byte, n uintptr) {
	var ptr uintptr
	ptr = uintptr(s)
	var i uintptr
	for i = 0; i < n; i++ {
		pByte := (*byte)(unsafe.Pointer(ptr + i))
		*pByte = c
	}
}
