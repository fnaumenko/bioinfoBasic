/**********************************************************
spline 2023 Fedor Naumenko (fedor.naumenko@gmail.com)
-------------------------
Last modified: 11/18/2023
-------------------------
Two-modes spline (smoothing curve) based in moving window
***********************************************************/
#pragma once

#include <vector>
#include <algorithm>	// sort
#include <memory>		// unique_ptr

enum eCurveType {
	SMOOTH = 0,		// sequential averaging by moving median and average
	SPIKED			// averaging by moving average only
};

// sliding spliner
template<typename T>	// T - type of smoothed values, assumed to be unsigned integer
class SSpliner {
	using splinesum_t = uint64_t;	// float is more versatile - for the future
	using slen_t = uint16_t;		// length type of moving window

public:
	// Returns 'silent zone length' - number of values pushed in after which the output is non-zero
	//	@param type: curve type
	//	@param base: half-length of moving window
	static slen_t SilentLength(eCurveType type, slen_t base) { return base * (2 + type) - 1; }

	// Constructor
	//	@param type: curve type
	//	@param base: half-length of moving window
	SSpliner(eCurveType type, slen_t base) :
		_baseLen(base), _curveType(type),
		_silentLen(SilentLength(type, base))
	{
		_ma.Init(base);
		if (_curveType) {
			_mm.reset(new MM<T>(base));
			_push = &MM<T>::Push;
		}
	}

	// Adds value and returns combined average
	//	@param val: input raw value
	float Push(T val) {
		_filledLen += (_filledLen <= _silentLen);
		return _ma.Push((_mm.get()->*_push)(val, _filledLen < _baseLen), _filledLen < _silentLen);
	}

	// Returns true X position
	chrlen CorrectX(chrlen x) const { return x - (_baseLen << _curveType); }

	// Returns 'silent zone length' - number of values pushed in after which the output is non-zero
	slen_t SilentLength() const { return _silentLen; }

private:
	// Moving window (Sliding subset)
	template<typename T>
	class MW : protected std::vector<T>
	{
	protected:
		MW() {}

		// Constructor by half of moving window length ("base")
		MW(slen_t base) { Init(base); }

		// Adds last value and pops the first one (QUEUE functionality)
		void PushVal(T val)
		{
			std::move(this->begin() + 1, this->end(), this->begin());
			*(this->end() - 1) = val;
		}

	public:
		// Initializes the instance to zeros
		//	@param base: half-length of moving window
		void Init(slen_t base) { this->insert(this->begin(), size_t(base) * 2 + 1, 0); }
	};

	// Simple Moving Average spliner
	// https://en.wikipedia.org/wiki/Moving_average
	template<typename T>
	class MA : public MW<T>
	{
		splinesum_t	_sum = 0;		// sum of adding values

	public:
		// Adds value and returns average
		//	@param val: input raw value
		//	@param zeroOutput: if true than return 0 (silent zone)
		float Push(T val, bool zeroOutput)
		{
			_sum += int64_t(val) - *this->begin();
			this->PushVal(val);
			return zeroOutput ? 0 : float(_sum) / this->size();
		}
	};

	// Simple Moving Median spliner
	template<typename T>
	class MM : protected MW<T>
	{
		std::vector<T> _ss;		// sorted moving window (sliding subset)

	public:
		// Constructor
		//	@param base: half-length of moving window
		MM(slen_t base) : MW<T>(base) { _ss.insert(_ss.begin(), this->size(), 0); }

		// Adds value and returns median
		//	@param val: input raw value
		//	@param zeroOutput: if true than return 0 (silent zone)
		T Push(T val, bool zeroOutput)
		{
			this->PushVal(val);
			if (zeroOutput)	return 0;
			std::copy(this->begin(), this->end(), _ss.begin());
			std::sort(_ss.begin(), _ss.end());
			return _ss[this->size() >> 1];		// mid-size
		}

		// Empty (stub) 'Push' method
		T PushStub(T val, bool) { return val; }
	};

	const eCurveType  _curveType;
	const slen_t	  _baseLen;			// half of moving window length
	const slen_t	  _silentLen;		// length of inner buffer that must be filled before the first valid output
	slen_t			  _filledLen = 0;	// number of filled values; limited by _silentLen
	MA<T>			  _ma;
	unique_ptr<MM<T>> _mm;
	T(MM<T>::* _push)(T, bool) = &MM<T>::PushStub;	// pointer to MM:Push function: real, or empty (stub) by default
};
