#include "RecoLocalCalo/EcalRecAlgos/interface/EcalUncalibRecHitMultiFitAlgo.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"

EcalUncalibRecHitMultiFitAlgo::EcalUncalibRecHitMultiFitAlgo() :
  _computeErrors(true),
  _doPrefit(false),
  _prefitMaxChiSq(1.0),
  _dynamicPedestals(false),
  _mitigateBadSamples(false),
  _selectiveBadSampleCriteria(false),
  _addPedestalUncertainty(0.),
  _simplifiedNoiseModelForGainSwitch(true),
  _gainSwitchUseMaxSample(false){

  _singlebx.resize(1);
  _singlebx << 0;

  _pulsefuncSingle.disableErrorCalculation();
  _pulsefuncSingle.setMaxIters(1);
  _pulsefuncSingle.setMaxIterWarnings(false);

}

/// compute rechits
EcalUncalibratedRecHit EcalUncalibRecHitMultiFitAlgo::makeRecHit(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const SampleMatrixGainArray &noisecors, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const BXVector &activeBX) {
  uint32_t flags = 0;

  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;

  double maxamplitude = -std::numeric_limits<double>::max();
  const unsigned int iSampleMax = 5;
  const unsigned int iFullPulseMax = 9;

  double pedval = 0.;

  SampleVector amplitudes;
  SampleGainVector gainsNoise;
  SampleGainVector gainsPedestal;
  SampleGainVector badSamples = SampleGainVector::Zero();
  bool hasSaturation = dataFrame.isSaturated();
  bool hasGainSwitch = hasSaturation || dataFrame.hasSwitchToGain6() || dataFrame.hasSwitchToGain1();

  //no dynamic pedestal in case of gain switch, since then the fit becomes too underconstrained
  bool dynamicPedestal = _dynamicPedestals && !hasGainSwitch;

  for(unsigned int iSample = 0; iSample < nsample; iSample++) {

    const EcalMGPASample &sample = dataFrame.sample(iSample);

    double amplitude = 0.;
    int gainId = sample.gainId();

    double pedestal = 0.;
    double gainratio = 1.;

    if (gainId==0 || gainId==3) {
      pedestal = aped->mean_x1;
      gainratio = aGain->gain6Over1()*aGain->gain12Over6();
      gainsNoise[iSample] = 2;
      gainsPedestal[iSample] = dynamicPedestal ? 2 : -1;  //-1 for static pedestal
    }
    else if (gainId==1) {
      pedestal = aped->mean_x12;
      gainratio = 1.;
      gainsNoise[iSample] = 0;
      gainsPedestal[iSample] = dynamicPedestal ? 0 : -1; //-1 for static pedestal
    }
    else if (gainId==2) {
      pedestal = aped->mean_x6;
      gainratio = aGain->gain12Over6();
      gainsNoise[iSample] = 1;
      gainsPedestal[iSample] = dynamicPedestal ? 1 : -1; //-1 for static pedestals
    }

    if (dynamicPedestal) {
      amplitude = (double)(sample.adc())*gainratio;
    }
    else {
      amplitude = ((double)(sample.adc()) - pedestal) * gainratio;
    }

    if (gainId == 0) {
       edm::LogError("EcalUncalibRecHitMultiFitAlgo")<< "Saturation encountered.  Multifit is not intended to be used for saturated channels.";
      //saturation
      if (dynamicPedestal) {
        amplitude = 4095.*gainratio;
      }
      else {
        amplitude = (4095. - pedestal) * gainratio;
      }
    }

    amplitudes[iSample] = amplitude;

    if (iSample==iSampleMax) {
      maxamplitude = amplitude;
      pedval = pedestal;
    }

  }

  double amplitude, amperr, chisq;
  bool status = false;

  //special handling for gain switch, where sample before maximum is potentially affected by slew rate limitation
  //optionally apply a stricter criteria, assuming slew rate limit is only reached in case where maximum sample has gain switched but previous sample has not
  //option 1: use simple max-sample algorithm
  if (hasGainSwitch && _gainSwitchUseMaxSample) {
    double maxpulseamplitude = maxamplitude / fullpulse[iFullPulseMax];
    EcalUncalibratedRecHit rh( dataFrame.id(), maxpulseamplitude, pedval, 0., 0., flags );
    rh.setAmplitudeError(0.);
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0) {
        rh.setOutOfTimeAmplitude(bx+5, 0.0);
      }
    }
    return rh;
  }

  //option2: A floating negative single-sample offset is added to the fit
  //such that the affected sample is treated only as a lower limit for the true amplitude
  bool mitigateBadSample = _mitigateBadSamples && hasGainSwitch && iSampleMax>0;
  mitigateBadSample &= (!_selectiveBadSampleCriteria || (gainsNoise.coeff(iSampleMax-1)!=gainsNoise.coeff(iSampleMax)) );
  if (mitigateBadSample) {
    badSamples[iSampleMax-1] = 1;
  }

  //compute noise covariance matrix, which depends on the sample gains
  SampleMatrix noisecov;
  if (hasGainSwitch) {
    std::array<double,3> pedrmss = {{aped->rms_x12, aped->rms_x6, aped->rms_x1}};
    std::array<double,3> gainratios = {{ 1., aGain->gain12Over6(), aGain->gain6Over1()*aGain->gain12Over6()}};
    if (_simplifiedNoiseModelForGainSwitch) {
      int gainidxmax = gainsNoise[iSampleMax];
      noisecov = gainratios[gainidxmax]*gainratios[gainidxmax]*pedrmss[gainidxmax]*pedrmss[gainidxmax]*noisecors[gainidxmax];
      if (!dynamicPedestal && _addPedestalUncertainty>0.) {
        //add fully correlated component to noise covariance to inflate pedestal uncertainty
        noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
      }
    }
    else {
      noisecov = SampleMatrix::Zero();
      for (unsigned int gainidx=0; gainidx<noisecors.size(); ++gainidx) {
        SampleGainVector mask = gainidx*SampleGainVector::Ones();
        SampleVector pedestal = (gainsNoise.array()==mask.array()).cast<SampleVector::value_type>();
        if (pedestal.maxCoeff()>0.) {
          //select out relevant components of each correlation matrix, and assume no correlation between samples with
          //different gain
          noisecov += gainratios[gainidx]*gainratios[gainidx]*pedrmss[gainidx]*pedrmss[gainidx]*pedestal.asDiagonal()*noisecors[gainidx]*pedestal.asDiagonal();
          if (!dynamicPedestal && _addPedestalUncertainty>0.) {
            //add fully correlated component to noise covariance to inflate pedestal uncertainty
            noisecov += gainratios[gainidx]*gainratios[gainidx]*_addPedestalUncertainty*_addPedestalUncertainty*pedestal.asDiagonal()*SampleMatrix::Ones()*pedestal.asDiagonal();
          }
        }
      }
    }
  }
  else {
    noisecov = aped->rms_x12*aped->rms_x12*noisecors[0];
    if (!dynamicPedestal && _addPedestalUncertainty>0.) {
      //add fully correlated component to noise covariance to inflate pedestal uncertainty
      noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
    }
  }

  //optimized one-pulse fit for hlt
  bool usePrefit = false;
  if (_doPrefit) {
    status = _pulsefuncSingle.DoFit(amplitudes,noisecov,_singlebx,fullpulse,fullpulsecov,gainsPedestal,badSamples);
    amplitude = status ? _pulsefuncSingle.X()[0] : 0.;
    amperr = status ? _pulsefuncSingle.Errors()[0] : 0.;
    chisq = _pulsefuncSingle.ChiSq();

    if (chisq < _prefitMaxChiSq) {
      usePrefit = true;
    }
  }

  if (!usePrefit) {

    if(!_computeErrors) _pulsefunc.disableErrorCalculation();
    status = _pulsefunc.DoFit(amplitudes,noisecov,activeBX,fullpulse,fullpulsecov,gainsPedestal,badSamples);
    chisq = _pulsefunc.ChiSq();


    if (!status) {
      edm::LogWarning("EcalUncalibRecHitMultiFitAlgo::makeRecHit") << "Failed Fit" << std::endl;
    }

    unsigned int ipulseintime = 0;
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      if (_pulsefunc.BXs().coeff(ipulse)==0) {
        ipulseintime = ipulse;
        break;
      }
    }

    amplitude = status ? _pulsefunc.X()[ipulseintime] : 0.;
    amperr = status ? _pulsefunc.Errors()[ipulseintime] : 0.;

  }

  double jitter = 0.;

  EcalUncalibratedRecHit rh( dataFrame.id(), amplitude , pedval, jitter, chisq, flags );
  rh.setAmplitudeError(amperr);

  if (!usePrefit) {
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0 && std::abs(bx)<100) {
        rh.setOutOfTimeAmplitude(bx+5, status ? _pulsefunc.X().coeff(ipulse) : 0.);
      }
      else if (bx==(100+gainsPedestal[iSampleMax])) {
        rh.setPedestal(status ? _pulsefunc.X().coeff(ipulse) : 0.);
      }
    }
  }

  return rh;
}

double EcalUncalibRecHitMultiFitAlgo::computeTime(const EcalDataFrame& dataFrame, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const SampleMatrixGainArray &noisecors, const FullSampleVector &fullpulse, const FullSampleMatrix &fullpulsecov, const BXVector &activeBX, const float startTime, const float stopTime, const float stepTime) {
  uint32_t flags = 0;

  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;

  double maxamplitude = -std::numeric_limits<double>::max();
  const unsigned int iSampleMax = 5;
  const unsigned int iFullPulseMax = 9;

  double pedval = 0.;

  SampleVector amplitudes;
  SampleGainVector gainsNoise;
  SampleGainVector gainsPedestal;
  SampleGainVector badSamples = SampleGainVector::Zero();
  bool hasSaturation = dataFrame.isSaturated();
  bool hasGainSwitch = hasSaturation || dataFrame.hasSwitchToGain6() || dataFrame.hasSwitchToGain1();

  //no dynamic pedestal in case of gain switch, since then the fit becomes too underconstrained
  bool dynamicPedestal = _dynamicPedestals && !hasGainSwitch;

  for(unsigned int iSample = 0; iSample < nsample; iSample++) {

    const EcalMGPASample &sample = dataFrame.sample(iSample);

    double amplitude = 0.;
    int gainId = sample.gainId();

    double pedestal = 0.;
    double gainratio = 1.;

    if (gainId==0 || gainId==3) {
      pedestal = aped->mean_x1;
      gainratio = aGain->gain6Over1()*aGain->gain12Over6();
      gainsNoise[iSample] = 2;
      gainsPedestal[iSample] = dynamicPedestal ? 2 : -1;  //-1 for static pedestal
    }
    else if (gainId==1) {
      pedestal = aped->mean_x12;
      gainratio = 1.;
      gainsNoise[iSample] = 0;
      gainsPedestal[iSample] = dynamicPedestal ? 0 : -1; //-1 for static pedestal
    }
    else if (gainId==2) {
      pedestal = aped->mean_x6;
      gainratio = aGain->gain12Over6();
      gainsNoise[iSample] = 1;
      gainsPedestal[iSample] = dynamicPedestal ? 1 : -1; //-1 for static pedestals
    }

    if (dynamicPedestal) {
      amplitude = (double)(sample.adc())*gainratio;
    }
    else {
      amplitude = ((double)(sample.adc()) - pedestal) * gainratio;
    }

    if (gainId == 0) {
       edm::LogError("EcalUncalibRecHitMultiFitAlgo")<< "Saturation encountered.  Multifit is not intended to be used for saturated channels.";
      //saturation
      if (dynamicPedestal) {
        amplitude = 4095.*gainratio;
      }
      else {
        amplitude = (4095. - pedestal) * gainratio;
      }
    }

    amplitudes[iSample] = amplitude;

    if (iSample==iSampleMax) {
      maxamplitude = amplitude;
      pedval = pedestal;
    }

  }

  double chisq;

  //special handling for gain switch, where sample before maximum is potentially affected by slew rate limitation
  //optionally apply a stricter criteria, assuming slew rate limit is only reached in case where maximum sample has gain switched but previous sample has not
  //option 1: use simple max-sample algorithm
  if (hasGainSwitch && _gainSwitchUseMaxSample) {
    double maxpulseamplitude = maxamplitude / fullpulse[iFullPulseMax];
    EcalUncalibratedRecHit rh( dataFrame.id(), maxpulseamplitude, pedval, 0., 0., flags );
    rh.setAmplitudeError(0.);
    for (unsigned int ipulse=0; ipulse<_pulsefunc.BXs().rows(); ++ipulse) {
      int bx = _pulsefunc.BXs().coeff(ipulse);
      if (bx!=0) {
        rh.setOutOfTimeAmplitude(bx+5, 0.0);
      }
    }
    return -100;
  }

  //option2: A floating negative single-sample offset is added to the fit
  //such that the affected sample is treated only as a lower limit for the true amplitude
  bool mitigateBadSample = _mitigateBadSamples && hasGainSwitch && iSampleMax>0;
  mitigateBadSample &= (!_selectiveBadSampleCriteria || (gainsNoise.coeff(iSampleMax-1)!=gainsNoise.coeff(iSampleMax)) );
  if (mitigateBadSample) {
    badSamples[iSampleMax-1] = 1;
  }

  //compute noise covariance matrix, which depends on the sample gains
  SampleMatrix noisecov;
  if (hasGainSwitch) {
    std::array<double,3> pedrmss = {{aped->rms_x12, aped->rms_x6, aped->rms_x1}};
    std::array<double,3> gainratios = {{ 1., aGain->gain12Over6(), aGain->gain6Over1()*aGain->gain12Over6()}};
    if (_simplifiedNoiseModelForGainSwitch) {
      int gainidxmax = gainsNoise[iSampleMax];
      noisecov = gainratios[gainidxmax]*gainratios[gainidxmax]*pedrmss[gainidxmax]*pedrmss[gainidxmax]*noisecors[gainidxmax];
      if (!dynamicPedestal && _addPedestalUncertainty>0.) {
        //add fully correlated component to noise covariance to inflate pedestal uncertainty
        noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
      }
    }
    else {
      noisecov = SampleMatrix::Zero();
      for (unsigned int gainidx=0; gainidx<noisecors.size(); ++gainidx) {
        SampleGainVector mask = gainidx*SampleGainVector::Ones();
        SampleVector pedestal = (gainsNoise.array()==mask.array()).cast<SampleVector::value_type>();
        if (pedestal.maxCoeff()>0.) {
          //select out relevant components of each correlation matrix, and assume no correlation between samples with
          //different gain
          noisecov += gainratios[gainidx]*gainratios[gainidx]*pedrmss[gainidx]*pedrmss[gainidx]*pedestal.asDiagonal()*noisecors[gainidx]*pedestal.asDiagonal();
          if (!dynamicPedestal && _addPedestalUncertainty>0.) {
            //add fully correlated component to noise covariance to inflate pedestal uncertainty
            noisecov += gainratios[gainidx]*gainratios[gainidx]*_addPedestalUncertainty*_addPedestalUncertainty*pedestal.asDiagonal()*SampleMatrix::Ones()*pedestal.asDiagonal();
          }
        }
      }
    }
  }
  else {
    noisecov = aped->rms_x12*aped->rms_x12*noisecors[0];
    if (!dynamicPedestal && _addPedestalUncertainty>0.) {
      //add fully correlated component to noise covariance to inflate pedestal uncertainty
      noisecov += _addPedestalUncertainty*_addPedestalUncertainty*SampleMatrix::Ones();
    }
  }



  if(!_computeErrors) _pulsefunc.disableErrorCalculation();
  _pulsefunc.DoFit(amplitudes,noisecov,activeBX,fullpulse,fullpulsecov,gainsPedestal,badSamples);
  chisq = _pulsefunc.ChiSq();
  // std::cout<<"NICKLOG PULSE: ";
  // for (int i=0;i<fullpulse.size(); ++i) {
  //   std::cout<<fullpulse[i]<<" ";
  // }
  // std::cout<<std::endl;
  double minChisq=chisq;
  double tMin = -100;

  for (double t = startTime; t < stopTime; t += stepTime) {
    _pulsefunc.DoFit(amplitudes,noisecov,activeBX,interpolate(fullpulse, t),fullpulsecov,gainsPedestal,badSamples);
    chisq = _pulsefunc.ChiSq();
    // std::cout<<"NICKLOG: shifted fit: shift "<< t <<" ns\tChiSq "<< chisq << std::endl;
    if (chisq < minChisq) {
      // std::cout<<"NICKLOG: "<< chisq <<" < "<< minChisq << std::endl;
      minChisq = chisq;
      tMin = t;
    }
    // std::cout<<"NICKLOG PULSE tshift "<< t <<" ns: ";
  //   auto ipulse = interpolate(fullpulse, t);
  //   for (int i=0;i<ipulse.size(); ++i) {
  //     std::cout<<ipulse[i]<<" ";
  //   }
  //   std::cout<<std::endl;
  }

  // std::cout<<"NICKLOG: selected fit: shift "<< tMin <<" ns\tChiSq "<< minChisq << "\t vs prev chisq "<< chisqPrev << std::endl;
  return tMin;
}

double EcalUncalibRecHitMultiFitAlgo::computeTimeCC(const EcalDataFrame& dataFrame, const std::vector<double> &amplitudes, const EcalPedestals::Item * aped, const EcalMGPAGainRatio * aGain, const FullSampleVector &fullpulse, const float startTime, const float stopTime, const float stepTime) {
  const unsigned int nsample = EcalDataFrame::MAXSAMPLES;

  double maxamplitude = -std::numeric_limits<double>::max();

  double pulsenorm = 0.;

  std::vector<double> pedSubSamples(nsample);
  for(unsigned int iSample = 0; iSample < nsample; iSample++) {

    const EcalMGPASample &sample = dataFrame.sample(iSample);

    double amplitude = 0.;
    int gainId = sample.gainId();

    double pedestal = 0.;
    double gainratio = 1.;

    if (gainId==0 || gainId==3) {
      pedestal = aped->mean_x1;
      gainratio = aGain->gain6Over1()*aGain->gain12Over6();
    }
    else if (gainId==1) {
      pedestal = aped->mean_x12;
      gainratio = 1.;
    }
    else if (gainId==2) {
      pedestal = aped->mean_x6;
      gainratio = aGain->gain12Over6();
    }

    amplitude = ((double)(sample.adc()) - pedestal) * gainratio;

    if (gainId == 0) {
      //saturation
      amplitude = (4095. - pedestal) * gainratio;
    }

    pedSubSamples.at(iSample) = amplitude;

    if (amplitude>maxamplitude) {
      maxamplitude = amplitude;
    }
    pulsenorm += fullpulse(iSample);
  }

  // std::vector<double>::const_iterator amplit;
  // for(amplit=amplitudes.begin(); amplit<amplitudes.end(); ++amplit) {
  //   int ipulse = std::distance(amplitudes.begin(),amplit);
  //   int bx = ipulse - 5;
  //   int firstsamplet = std::max(0,bx + 3);
  //   int offset = 7-3-bx;
  //
  //   TVectorD pulse;
  //   pulse.ResizeTo(nsample);
  //   for (unsigned int isample = firstsamplet; isample<nsample; ++isample) {
  //     pulse(isample) = fullpulse(isample+offset);
  //     pedSubSamples(isample) = std::max(0., pedSubSamples(isample) - amplitudes[ipulse]*pulse(isample)/pulsenorm);
  //   }
  // }



  // std::cout<<"NICKLOG input: ";
  // for (unsigned int i=0;i<pedSubSamples.size(); ++i) {
  //   std::cout<<pedSubSamples[i]<<" ";
  // }
  // std::cout<<std::endl;

  double minDist=100000;
  double tMin = -100;

  //find shift:
  int shift=0;
  for (int s=0; s<10; ++s) {
    double dist = .0;
    for (unsigned int i=0; i<pedSubSamples.size(); ++i){
      if (i+s>fullpulse.size()) break;
      dist += std::pow(fullpulse[i+s]-pedSubSamples[i],2);
    }
    if (dist < minDist) {
      minDist = dist;
      shift = s;
    }
  }
  // std::cout<<"NICKLOG: shift "<< shift << std::endl;

  for (double t = startTime; t < stopTime; t += stepTime) {
    auto interpolated = interpolate(fullpulse, t);
    double dist = .0;

    for (unsigned int i=0; i<pedSubSamples.size(); ++i){
      if (i+shift>interpolated.size()) break;
      dist += std::pow(interpolated[i+shift]-pedSubSamples[i],2);
    }
    // std::cout<<"NICKLOG: t: "<< t <<" ns\tdist "<< dist << std::endl;
    if (dist < minDist) {
      // std::cout<<"NICKLOG: "<< chisq <<" < "<< minChisq << std::endl;
      minDist = dist;
      tMin = t;
    }
    // std::cout<<"NICKLOG PULSE tshift "<< t <<" ns: ";
  //   auto ipulse = interpolate(fullpulse, t);
  //   for (int i=0;i<ipulse.size(); ++i) {
  //     std::cout<<ipulse[i]<<" ";
  //   }
  //   std::cout<<std::endl;
  }

  // std::cout<<"NICKLOG: selected fit: shift "<< tMin <<" ns\tChiSq "<< minChisq << "\t vs prev chisq "<< chisqPrev << std::endl;
  return tMin + (shift-5)*25.;
}

FullSampleVector EcalUncalibRecHitMultiFitAlgo::interpolate(const FullSampleVector& fullpulse, const float t){
  FullSampleVector interpPulse;
  interpPulse[0] = fullpulse[0];
  interpPulse[interpPulse.size()-1] = fullpulse[interpPulse.size()-1];
  if (t==0) {
    for (int i=1; i<interpPulse.size()-1; ++i)
      interpPulse[i] = fullpulse[i];
  }
  else if (t>0) {
    for (int i=1; i<interpPulse.size()-1; ++i){
      interpPulse[i] = fullpulse[i] + (t/25.)*(fullpulse[i+1]-fullpulse[i]);
    }
  }
  else {
    for (int i=1; i<interpPulse.size()-1; ++i){
      interpPulse[i] = fullpulse[i-1] + (t/25.+1)*(fullpulse[i]-fullpulse[i-1]);
    }
  }
  return interpPulse;
}
