
package main
import (
    "image/color"
    "image/png"
    "os"
    "image"
    "fmt"
    "math"
)

// image output stuff

func newImg(width, height int) (*image.RGBA, func(int, int, int, int, int)) {
    img := image.NewRGBA(image.Rectangle{image.Point{0, 0}, image.Point{width, height}})
    return img, func(x, y, r, g, b int) {
        // img.Set(x, y, color.RGBA{uint8(rad % 256), uint8(grn % 256), uint8(blu % 256), 0xff})
        img.Set(x, y, color.RGBA{uint8(r % 256), uint8(g % 256), uint8(b % 256), 0xff})
    }
}

func dumpImg(img *image.RGBA) {
    f, err := os.Create(fileName())
    if err != nil {panic(err)}
    png.Encode(f, img) // encode as png
}

func fileName() string {
    for i := 0;; i++ {
        name := fmt.Sprintf("/root/0Git/fracGen/images/image%v.png", i)
        _, err := os.Stat(name)
        if err != nil {return name}
    }
}

// math

func mapRange(fs, fe, ts, te float64) func(float64) float64 {
    scale := ((te-ts)/(fe-fs))
    toff := (te+ts)/2
    foff := (fe+fs)/2
    return func(num float64) float64 {
        // return (num + (ts+te)/2 - (fs+fe)/2)*((te-ts)/(fe-fs))
        return (num-foff)*scale + toff
    }
}

func chopRange(s, e float64) func(float64) float64 {
    return func(num float64) float64 {
        if num < s {return s} else if num > e {return e} else {return num}
    }
}

// returns a square range around x, y with zoom x and y ranges as powers of 2
func xyrange(pow, x, y float64) (float64, float64, float64, float64) {
    nudge := math.Exp2(pow-1)
    return x-nudge, x+nudge, y+nudge, y-nudge // y+nudge first cuz the y is flipped in the computer things or something
}
