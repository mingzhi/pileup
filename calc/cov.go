package calc

// Covariance contains cov structure
type Covariance struct {
	XY, X, Y float64
	N        int
}

// NewCovariance create a Covariance
func NewCovariance() *Covariance {
	return &Covariance{}
}

// Increment add data to the calculator
func (c *Covariance) Increment(x, y float64) {
	c.XY += x * y
	c.X += x
	c.Y += y
	c.N++
}

// Append merges another covariance.
func (c *Covariance) Append(c1 *Covariance) {
	c.XY += c1.XY
	c.X += c1.X
	c.Y += c1.Y
	c.N += c1.N
}

// GetResult returns the result.
func (c *Covariance) GetResult() float64 {
	var v float64
	v = c.XY/float64(c.N) - (c.X/float64(c.N))*(c.Y/float64(c.N))
	return v
}

// GetN returns N.
func (c *Covariance) GetN() int {
	return c.N
}

// GetMeanX return mean of X.
func (c *Covariance) GetMeanX() float64 {
	return c.X / float64(c.N)
}

// GetMeanY returns mean of Y.
func (c *Covariance) GetMeanY() float64 {
	return c.Y / float64(c.N)
}