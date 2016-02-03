package calc

import (
	"github.com/mingzhi/gomath/stat/correlation"
)

// Calculator contains individual calculators.
type Calculator struct {
	MaxL int
	Cs   *MeanVariances
	Cr   *Covariances
	Ct   *Covariances
}

// New returns a new Calculator
func New(maxl int) *Calculator {
	c := Calculator{}
	c.MaxL = maxl
	c.Cs = NewMeanVariances(maxl)
	c.Cr = NewCovariances(maxl)
	c.Ct = NewCovariances(maxl)
	return &c
}

// Increment add x and y arrays, which separate at distance l.
func (c *Calculator) Increment(xArr, yArr []float64, l int) {
	if l < c.MaxL {
        // calculate covariance of x and y.
		cov := correlation.NewBivariateCovariance(false)
		for i := range xArr {
			x, y := xArr[i], yArr[i]
			cov.Increment(x, y)
		}
        
		c.Cs.Increment(l, cov.GetResult())
		c.Ct.AppendAt(l, cov)
		c.Cr.Increment(l, cov.MeanX(), cov.MeanY())

	}
}

// Append appends a calculator to another
func (c *Calculator) Append(c1 *Calculator) {
	c.Cs.Append(c1.Cs)
	c.Cr.Append(c1.Cr)
	c.Ct.Append(c1.Ct)
}