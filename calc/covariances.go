package calc


// Covariances contains an array of Coveriance.
type Covariances struct {
	corrs []*Covariance
}

// NewCovariances create a new Covariances
func NewCovariances(maxl int) *Covariances {
	cc := Covariances{}
	for i := 0; i < maxl; i++ {
		bc := NewCovariance()
		cc.corrs = append(cc.corrs, bc)
	}
	return &cc
}

// Increment add data (x, y) to the l Covariance.
func (c *Covariances) Increment(l int, x, y float64) {
	c.corrs[l].Increment(x, y)
}

// Append append a Covariance to the l Covariance.
func (c *Covariances) Append(c1 *Covariances) {
	for i := 0; i < len(c.corrs); i++ {
		c.corrs[i].Append(c1.corrs[i])
	}
}

func (c *Covariances) AppendAt(i int, c1 *Covariance) {
	c.corrs[i].Append(c1)
}

// GetResult returns the result.
func (c *Covariances) GetResult(l int) float64 {
	return c.corrs[l].GetResult()
}

// GetN returns the number of data points.
func (c *Covariances) GetN(l int) int {
	return c.corrs[l].GetN()
}

// GetMeanX returns x_bar.
func (c *Covariances) GetMeanX(l int) float64 {
	return c.corrs[l].GetMeanX()
}

// GetMeanY returns y_bar.
func (c *Covariances) GetMeanY(l int) float64 {
	return c.corrs[l].GetMeanY()
}
