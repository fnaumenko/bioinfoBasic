#pragma once

#include <vector>
#include <algorithm>	// sort
#include <memory>		// unique_ptr

// sliding spliner
class SSpliner {
public: 
	using slen_t = uint16_t;	// length of moving window

private:
	using spline_t = uint32_t;
	using splinesum_t = uint64_t;	// float;

	// Moving window (Sliding subset)
	class MW : protected std::vector<spline_t>
	{
	protected:
		MW() {}

		// Constructor by half of moving window length ("base")
		MW(slen_t base) { Init(base); }

		// Add last value and pop the first one (QUEUE functionality)
		void PushVal(spline_t val)
		{
			std::move(begin() + 1, end(), begin());
			*(end() - 1) = val;
		}

	public:
		void Init(slen_t base) { insert(begin(), size_t(base) * 2 + 1, 0); }
	};

	// Simple Moving Average splicer
	// https://en.wikipedia.org/wiki/Moving_average
	class MA : public MW
	{
		splinesum_t	_sum = 0;		// sum of adding values

	public:
		// Add value and return average
		//	@param val: input raw value
		//	@param zeroOutput: if true than return 0 (silent zone)
		float Push(spline_t val, bool zeroOutput)
		{
			_sum += val - *begin();
			PushVal(val);
			return zeroOutput ? 0 : float(_sum) / size();
		}
	};

	// Simple Moving Median splicer
	class MM : protected MW
	{
		std::vector<spline_t> _ss;		// sorted moving window (sliding subset)

	public:
		// Constructor
		//	@param base: half-length of moving window
		MM(slen_t base) : MW(base) { _ss.insert(_ss.begin(), size(), 0); }

		// Add value and return median
		//	@param val: input raw value
		//	@param zeroOutput: if true than return 0 (silent zone)
		spline_t Push(spline_t val, bool zeroOutput)
		{
			PushVal(val);
			if (zeroOutput)	return 0;
			std::copy(begin(), end(), _ss.begin());
			std::sort(_ss.begin(), _ss.end());
			return _ss[size() >> 1];		// mid-size
		}

		// Empty (stub) 'Push' method
		spline_t PushStub(spline_t val, bool) { return val; }
	};

public:
	enum eCurveType { SMOOTH = 0, SPIKED };
	//const static slen_t ReadSplineBase = 10;	//half-length of moving window for reads spline

	static slen_t SilentLength(eCurveType type, slen_t base)// = ReadSplineBase)
	{
		return base * 2 + base * type - 1;
	}

	// Constructor
	//	@param type: curve type
	//	@param base: half-length of moving window
	SSpliner(eCurveType type, slen_t base/* = ReadSplineBase*/) :
		_base(base), _curveType(type),
		_silentLen(SilentLength(type, base))
	{
		_ma.Init(base);
		if (_curveType) {
			_mm.reset(new MM(base));
			_push = &MM::Push;
		}
	}

	// Add value and return combined average
	//	@param val: input raw value
	float Push(spline_t val) {
		_valCnt += (_valCnt <= _silentLen);
		//if (_valCnt <= _silentLen)	_valCnt++;
		return _ma.Push((_mm.get()->*_push)(val, _valCnt <= _base), _valCnt < _silentLen);
	}

	// Returns true X position
	chrlen CorrectX(chrlen x) const { return x - (_base << _curveType); }

	// Returns length of inner buffer that must be filled before the first valid output 
	spline_t SilentLength() const { return _silentLen; }

private:
	const eCurveType	_curveType;
	const slen_t	_base;			// half of moving window length
	const slen_t	_silentLen;		// length of inner buffer that must be filled before the first valid output
	slen_t			_valCnt = 0;	// value counter, limited by _silentLen
	MA				_ma;
	unique_ptr<MM>	_mm;
	spline_t(MM::* _push)(spline_t, bool) = &MM::PushStub;	// pointer to MM:Push function: real or empty (stub)
};
