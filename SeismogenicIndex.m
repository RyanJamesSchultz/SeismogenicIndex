function [Veq, SIGMA,sigma,sigma_err, R2, Vfit,Nfit, Vs]=SeismogenicIndex(Tv, Vc, Tm, M, b, Mc, Vstart,Vend)
  % SeismogenicIndex(Tv, Vc, Tm, M, b, Mc, Vstart,Vend)
  %
  % A function that calculates the Seismogenic Index as described in 
  % Shapiro et al. (2007; 2010).
  %
  % Input data:
  % Tv   - Cumulative volume injected time axis (same units as Tm).
  % Vc   - Cumulative volume injected (m³).
  % Tm   - Earthquake time axis (same units as Tv).
  % M    - Earthquake magnitudes (same units as Mc).
  % b    - Seismic b-value.
  % Mc   - Magnitude cut-off (same units as M).
  %
  % Output data:
  % Veq        - Cumulative injected volume at the time of each earthquake (m³).
  % SIGMA      - Seismogenic Index as a function of time/volume.
  % sigma      - Best guess at a single Seismogenic Index value.
  % sigma_err  - Variability in single Seismogenic Index value.
  % R2         - Goodness-of-fit statistic.
  % Vfit       - Cumulative volume fit values (m³), for plotting.
  % Nfit       - Cumulative earthquake count fit values, for plotting.
  % Vs         - Cumulative injected volume (m³) at time of first truncated earthquake.
  %
  % References:
  %   Shapiro, Dinske, & Kummerow (2007) Geophys. Res. Lett. 34(22), doi: 10.1029/2007GL031615.
  %   Shapiro, Dinske, Langenbruch, & Wenzel (2010) TLE. 29(3),304-309, doi: 10.1190/1.3353727.
  %
  %
  % This program is free software: you can redistribute it and/or modify
  % it under the terms of the GNU General Public License as published
  % by the Free Software Foundation, either version 3 of the License, or
  % any later version.
  % This program is distributed in the hope that it will be useful,
  % but WITHOUT ANY WARRANTY; without even the implied warranty of
  % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  % GNU General Public License for more details:
  % http://www.gnu.org/licenses/
  
  
  % Check dimensions and sizes of input vectors.
  if( (length(Tv)~=length(Vc))||(length(Tm)~=length(M))   )
      fprintf('improper input lengths.\n');
      Veq=0; SIGMA=0;
      return;
  end;
  if(~isrow(Tv))
      Tv=Tv';
  end;
  if(~isrow(Vc))
      Vc=Vc';
  end;
  if(~isrow(Tm))
      Tm=Tm';
  end;
  if(~isrow(M))
      M=M';
  end;
  
  % Sort earthquakes by time.
  [Tm, t_ind]=sort(Tm);
  M=M(t_ind);
  
  % Keep only earthquakes above given threshold.
  Tm=Tm(M>=Mc);
  M=M(M>=Mc);
  
  % Check to see there are still earthquakes to use.
  if(isempty(M))
      fprintf('No earthquakes above given magnitude threshold.\n');
      Veq=0; SIGMA=0;
      return;
  end;
  
  % Interpolate the volume of injectant at each earthquake time.
  Veq=interp1(Tv, Vc, Tm, 'linear');
  
  % Keep only earthquakes happening before user defined shut-in volume.
  if(Vend~=0)
      Tm=Tm(Veq<=Vend);  M=M(Veq<=Vend);  Veq=Veq(Veq<=Vend);
  end;
  
  % Keep only earthquakes happening after injection interval or user defined start volume.
  Vs=Veq(1);
  if(Vstart==0)
      Veq=Veq-Vs;
      Tm=Tm(Veq>0);   M=M(Veq>0);  Veq=Veq(Veq>0);
  else
      Tm=Tm(Veq>=Vstart);  M=M(Veq>=Vstart);  Veq=Veq(Veq>=Vstart);
      Veq=Veq-Vstart;
  end;
  
  % Keep only earthquakes that have interpolated injection volumes.
  Tm=Tm(~isnan(Veq));   M=M(~isnan(Veq));   Veq=Veq(~isnan(Veq));
  
  %check to see there are still earthquakes to use.
  if(isempty(M))
      fprintf('No earthquakes during injection interval.\n');
      Veq=0; SIGMA=0;
      return;
  end;
  
  % Calculate seismogenic index vector.
  Nm=1:length(M);
  Y=log10(Nm);
  X=log10(Veq);
  SIGMA=Y-X+b*Mc;
  
  % Best guess at single seismogenic index and its error.
  w=(Nm).^2; w=w/sum(w);
  sigma=sum(SIGMA.*w);
  sigma_err=sqrt(sum((SIGMA.^2).*w) - sigma^2);
  
  % Determine goodness-of-fit (in linear space).
  Y=Nm;
  Yfit=10.^(X-b*Mc+sigma);
  Ybar=mean(Y);
  SStot=sum(w.*((Y-Ybar).^2));
  SSres=sum(w.*((Y-Yfit).^2));
  R2=1-(SSres/SStot);
  
  % Fit to line values.
  Vfit=linspace(Veq(1), Veq(end), length(M));
  Nfit=10.^(log10(Vfit) - b*Mc + sigma );
  
return;


