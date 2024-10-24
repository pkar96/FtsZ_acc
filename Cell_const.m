classdef Cell_const < matlab.mixin.Copyable
   properties
      t
      v
      s
      n
      nm
      nthresh
      rate
      oris
      vNextInit
      vNextDivs
      tNextDivs
      tNextConst_in
      vNextConst
      vOfInits
      tOfInits
      oOfInits
      vOfConst
      vb
      vi
      vd
      tLastDiv
      tLastInit
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
      tot_cell
      rec_data
   end
   methods
       function obj = Cell_const(v)
         if nargin == 0
            obj.t = 0;
            obj.v = 1;
            obj.s = 1;
            obj.n = 50;
            obj.nm= 0.2*50;
            obj.nthresh = 360;
            obj.vb = 1;
            obj.vi = 1;
            obj.vd = 1;
            obj.rate = 0;
            obj.vNextInit = 0;
            obj.vNextDivs = nan;
            obj.tNextDivs = nan;
            obj.tNextConst_in= [];
            obj.vNextConst = nan;
            obj.oris = 2;
            obj.tLastDiv = 0;
            obj.tLastInit = 0;
            obj.vOfInits = [1];
            obj.tOfInits = [0];
            obj.oOfInits = [2];
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
            obj.tot_cell = 0;
            obj.rec_data=0;
         else
            obj.t = 0;
            obj.v = v;
            obj.s = 1;
            obj.n = v*50;
            obj.nm= 0.2*v*50;
            obj.nthresh = 360;
            obj.vb = v;
            obj.vi = v;
            obj.vd = v;
            obj.rate = 0;
            obj.vNextInit = 0;
            obj.vNextDivs = nan;
            obj.tNextDivs = nan;
            obj.tNextConst_in= [];
            obj.vNextConst = nan;
            obj.oris = 2;
            obj.tLastDiv = 0;
            obj.tLastInit = 0;
            obj.vOfInits = [v];
            obj.tOfInits = [0];
            obj.oOfInits = [2];
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
            obj.tot_cell = 0;
            obj.rec_data=0;
         end
      end
   end
end