package calc

import (
	"github.com/mingzhi/gomath/stat/desc/meanvar"
)


// MeanVariances is an array of MeanVar.
type MeanVariances struct {
	meanvars []*meanvar.MeanVar
}

// NewMeanVariances return a new MeanVariances of length maxl.
func NewMeanVariances(size int) *MeanVariances {
	mcc := MeanVariances{}
	for i := 0; i < size; i++ {
		mv := meanvar.New()
		mcc.meanvars = append(mcc.meanvars, mv)
	}
	return &mcc
}

// Increment add a data point to the lst MeanVar.
func (m *MeanVariances) Increment(l int, v float64) {
	m.meanvars[l].Increment(v)
}

// GetMean returns the mean of the lst MeanVar.
func (m *MeanVariances) GetMean(l int) float64 {
	return m.meanvars[l].Mean.GetResult()
}

// GetVar returns the variance of the lst MeanVar.
func (m *MeanVariances) GetVar(l int) float64 {
	return m.meanvars[l].Var.GetResult()
}

// GetN returns the size of the data points of the lst MeanVar.
func (m *MeanVariances) GetN(l int) int {
	return m.meanvars[l].Mean.GetN()
}

func (m *MeanVariances) Size() int {
	return len(m.meanvars)
}

// Appends append a MeanVariances to the other
func (m *MeanVariances) Append(m1 *MeanVariances) {
	for i := 0; i < len(m.meanvars); i++ {
		m.meanvars[i].Mean.Append(m1.meanvars[i].Mean)
		m.meanvars[i].Var.Append(m1.meanvars[i].Var)
	}
}
