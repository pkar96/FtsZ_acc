classdef Cell_const_div < matlab.mixin.Copyable
   properties
      t
      v
      s
      n
      nm
      nthresh
      s_1
      n_1
      nm_1
      nthresh_1
      rate
      vNextDivs
      tNextDivs
      vNextConst
      vb
      vd
      tLastDiv
      vOfConst
      td_cell
      lb_cell
      ld_cell
      tc_cell
      lc_cell
      c_tc_cell
      cm_tc_cell
      rate_cell
      cm_cell
      c_cell
      rep_cell
      ftsz_cell
      tot_cell
      rec_data
   end
   methods
       function obj = Cell_const_div(v)
         if nargin == 0
            obj.t = 0;
            obj.v = 1;
            obj.s = 1;
            obj.n = 50;
            obj.nm= 0.2*50;
            obj.nthresh = 360;
            obj.s_1 = 1;
            obj.n_1 = 50;
            obj.nm_1= 0.2*50;
            obj.nthresh_1 = 360;
            obj.vb = 1;
            obj.vd = 1;
            obj.rate = 0;
            obj.vNextDivs = nan;
            obj.tNextDivs = nan;
            obj.vNextConst = nan;
            obj.tLastDiv = 0;
            obj.vOfConst = [];
            obj.td_cell=NaN;
            obj.lb_cell= NaN;
            obj.ld_cell= NaN;
            obj.rate_cell= NaN;
            obj.cm_cell = NaN;
            obj.c_cell = NaN;
            obj.tc_cell= NaN;
            obj.lc_cell= NaN;
            obj.c_tc_cell=NaN;
            obj.cm_tc_cell=NaN;
            obj.rep_cell = 0;
            obj.ftsz_cell = 0;
            obj.tot_cell = 0;
            obj.rec_data=0;
         else
            obj.t = 0;
            obj.v = v;
            obj.s = 1;
            obj.n = v*50;
            obj.nm= 0.2*v*50;
            obj.nthresh = 360;
            obj.s_1 = 1;
            obj.n_1 = v*50;
            obj.nm_1= 0.2*v*50;
            obj.nthresh_1 = 360;
            obj.vb = v;
            obj.vd = v;
            obj.rate = 0;
            obj.vNextDivs = nan;
            obj.tNextDivs = nan;
            obj.vNextConst = nan;
            obj.tLastDiv = 0;
            obj.vOfConst = [];
            obj.td_cell= NaN;
            obj.lb_cell= NaN;
            obj.ld_cell= NaN;
            obj.rate_cell= NaN;
            obj.cm_cell = NaN;
            obj.c_cell = NaN;
            obj.tc_cell= NaN;
            obj.lc_cell= NaN;
            obj.c_tc_cell=NaN;
            obj.cm_tc_cell=NaN;
            obj.rep_cell = 0;
            obj.ftsz_cell = 0;
            obj.tot_cell = 0;
            obj.rec_data=0;
         end
      end
   end
end