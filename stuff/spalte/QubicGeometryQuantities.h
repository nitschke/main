/*
 * QubicGeometryQuantities.h
 *
 *  Created on: 21.05.2014
 *      Author: sebastian
 */

#ifndef QUBICGEOMETRYQUANTITIES_H_
#define QUBICGEOMETRYQUANTITIES_H_

class QubicGeometryQuantities
{
public:

	QubicGeometryQuantities() {}

	static inline double getN0(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double y6 = y[6] ;
		double y7 = y[7] ;
		double y8 = y[8] ;
		double y9 = y[9] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t3 = x1*x1;
		double t4 = x0*x0;
		double t5 = x2*x2;
		double t8 = t4*y0*3.0;
		double t9 = t3*y5*3.0;
		double t10 = t5*y7*3.0;
		double t11 = x0*x1*y3*6.0;
		double t12 = x0*x2*y4*6.0;
		double t13 = x1*x2*y6*6.0;
		double t2 = t8+t9+t10+t11+t12+t13;
		double t6 = t3*y1*3.0+t4*y3*3.0+t5*y9*3.0+x0*x1*y5*6.0+x0*x2*y6*6.0+x1*x2*y8*6.0;
		double t7 = t5*y2*3.0+t4*y4*3.0+t3*y8*3.0+x0*x1*y6*6.0+x0*x2*y7*6.0+x1*x2*y9*6.0;
		return t2*1.0/sqrt((t2*t2)*(1.0/9.0)+(t6*t6)*(1.0/9.0)+(t7*t7)*(1.0/9.0))*(1.0/3.0);
	}

	static inline double getN1(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double y6 = y[6] ;
		double y7 = y[7] ;
		double y8 = y[8] ;
		double y9 = y[9] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t16 = x1*x1;
		double t17 = x0*x0;
		double t18 = x2*x2;
		double t15 = t17*y0*3.0+t16*y5*3.0+t18*y7*3.0+x0*x1*y3*6.0+x0*x2*y4*6.0+x1*x2*y6*6.0;
		double t21 = t16*y1*3.0;
		double t22 = t17*y3*3.0;
		double t23 = t18*y9*3.0;
		double t24 = x0*x1*y5*6.0;
		double t25 = x0*x2*y6*6.0;
		double t26 = x1*x2*y8*6.0;
		double t19 = t21+t22+t23+t24+t25+t26;
		double t20 = t18*y2*3.0+t17*y4*3.0+t16*y8*3.0+x0*x1*y6*6.0+x0*x2*y7*6.0+x1*x2*y9*6.0;
		return t19*1.0/sqrt((t15*t15)*(1.0/9.0)+(t19*t19)*(1.0/9.0)+(t20*t20)*(1.0/9.0))*(1.0/3.0);
	}

	static inline double getN2(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double y6 = y[6] ;
		double y7 = y[7] ;
		double y8 = y[8] ;
		double y9 = y[9] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t29 = x1*x1;
		double t30 = x0*x0;
		double t31 = x2*x2;
		double t28 = t30*y0*3.0+t29*y5*3.0+t31*y7*3.0+x0*x1*y3*6.0+x0*x2*y4*6.0+x1*x2*y6*6.0;
		double t32 = t29*y1*3.0+t30*y3*3.0+t31*y9*3.0+x0*x1*y5*6.0+x0*x2*y6*6.0+x1*x2*y8*6.0;
		double t34 = t30*y4*3.0;
		double t35 = t31*y2*3.0;
		double t36 = t29*y8*3.0;
		double t37 = x0*x1*y6*6.0;
		double t38 = x0*x2*y7*6.0;
		double t39 = x1*x2*y9*6.0;
		double t33 = t34+t35+t36+t37+t38+t39;
		return t33*1.0/sqrt((t28*t28)*(1.0/9.0)+(t32*t32)*(1.0/9.0)+(t33*t33)*(1.0/9.0))*(1.0/3.0);
	}

	static inline double getLaplaceN0(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double y6 = y[6] ;
		double y7 = y[7] ;
		double y8 = y[8] ;
		double y9 = y[9] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t42 = x1*x1;
		double t43 = x0*x0;
		double t44 = x2*x2;
		double t47 = t43*y0*3.0;
		double t48 = t42*y5*3.0;
		double t49 = t44*y7*3.0;
		double t50 = x0*x1*y3*6.0;
		double t51 = x0*x2*y4*6.0;
		double t52 = x1*x2*y6*6.0;
		double t41 = t47+t48+t49+t50+t51+t52;
		double t55 = t42*y1*3.0;
		double t56 = t43*y3*3.0;
		double t57 = t44*y9*3.0;
		double t58 = x0*x1*y5*6.0;
		double t59 = x0*x2*y6*6.0;
		double t60 = x1*x2*y8*6.0;
		double t45 = t55+t56+t57+t58+t59+t60;
		double t63 = t43*y4*3.0;
		double t64 = t44*y2*3.0;
		double t65 = t42*y8*3.0;
		double t66 = x0*x1*y6*6.0;
		double t67 = x0*x2*y7*6.0;
		double t68 = x1*x2*y9*6.0;
		double t46 = t63+t64+t65+t66+t67+t68;
		double t53 = t41*t41;
		double t54 = t53*(1.0/9.0);
		double t61 = t45*t45;
		double t62 = t61*(1.0/9.0);
		double t69 = t46*t46;
		double t70 = t69*(1.0/9.0);
		double t71 = t54+t62+t70;
		double t72 = 1.0/t71;
		double t73 = t53*t72*(1.0/9.0);
		double t74 = t73-1.0;
		double t76 = x0*y0*6.0;
		double t77 = x1*y3*6.0;
		double t78 = x2*y4*6.0;
		double t79 = t76+t77+t78;
		double t80 = t41*t79*(2.0/9.0);
		double t81 = x0*y3*6.0;
		double t82 = x1*y5*6.0;
		double t83 = x2*y6*6.0;
		double t84 = t81+t82+t83;
		double t85 = t45*t84*(2.0/9.0);
		double t86 = x0*y4*6.0;
		double t87 = x1*y6*6.0;
		double t88 = x2*y7*6.0;
		double t89 = t86+t87+t88;
		double t90 = t46*t89*(2.0/9.0);
		double t75 = t80+t85+t90;
		double t91 = 1.0/pow(t71,3.0/2.0);
		double t92 = 1.0/sqrt(t71);
		double t93 = x1*y1*6.0;
		double t94 = x0*y5*6.0;
		double t95 = x2*y8*6.0;
		double t96 = t93+t94+t95;
		double t97 = x0*y6*6.0;
		double t98 = x1*y8*6.0;
		double t99 = x2*y9*6.0;
		double t100 = t97+t98+t99;
		double t101 = t41*t84*(2.0/9.0);
		double t102 = t45*t96*(2.0/9.0);
		double t103 = t46*t100*(2.0/9.0);
		double t104 = t101+t102+t103;
		double t105 = 1.0/pow(t71,5.0/2.0);
		double t106 = x2*y2*6.0;
		double t107 = x0*y7*6.0;
		double t108 = x1*y9*6.0;
		double t109 = t106+t107+t108;
		double t110 = t41*t89*(2.0/9.0);
		double t111 = t46*t109*(2.0/9.0);
		double t112 = t45*t100*(2.0/9.0);
		double t113 = t110+t111+t112;
		double t114 = t84*t92*(1.0/3.0);
		double t118 = t41*t91*t104*(1.0/6.0);
		double t115 = t114-t118;
		double t116 = t89*t92*(1.0/3.0);
		double t120 = t41*t91*t113*(1.0/6.0);
		double t117 = t116-t120;
		double t119 = 1.0/(t71*t71);
		double t121 = t61*t72*(1.0/9.0);
		double t122 = t121-1.0;
		double t123 = t84*t84;
		double t124 = t123*(2.0/9.0);
		double t125 = t79*t84*(2.0/9.0);
		double t126 = t84*t96*(2.0/9.0);
		double t127 = t89*t100*(2.0/9.0);
		double t128 = t41*y3*(4.0/3.0);
		double t129 = t45*y5*(4.0/3.0);
		double t130 = t46*y6*(4.0/3.0);
		double t131 = t125+t126+t127+t128+t129+t130;
		double t132 = t41*t91*t131*(1.0/6.0);
		double t133 = t75*t84*t91*(1.0/6.0);
		double t134 = t79*t91*t104*(1.0/6.0);
		double t167 = t92*y3*2.0;
		double t168 = t41*t75*t104*t105*(1.0/4.0);
		double t135 = t132+t133+t134-t167-t168;
		double t136 = t79*t92*(1.0/3.0);
		double t138 = t41*t75*t91*(1.0/6.0);
		double t137 = t136-t138;
		double t139 = t69*t72*(1.0/9.0);
		double t140 = t139-1.0;
		double t141 = t89*t89;
		double t142 = t141*(2.0/9.0);
		double t143 = t100*t100;
		double t144 = t143*(2.0/9.0);
		double t145 = t79*t89*(2.0/9.0);
		double t146 = t89*t109*(2.0/9.0);
		double t147 = t84*t100*(2.0/9.0);
		double t148 = t41*y4*(4.0/3.0);
		double t149 = t45*y6*(4.0/3.0);
		double t150 = t46*y7*(4.0/3.0);
		double t151 = t145+t146+t147+t148+t149+t150;
		double t152 = t41*t91*t151*(1.0/6.0);
		double t153 = t75*t89*t91*(1.0/6.0);
		double t154 = t79*t91*t113*(1.0/6.0);
		double t179 = t92*y4*2.0;
		double t180 = t41*t75*t105*t113*(1.0/4.0);
		double t155 = t152+t153+t154-t179-t180;
		double t156 = t84*t89*(2.0/9.0);
		double t157 = t96*t100*(2.0/9.0);
		double t158 = t100*t109*(2.0/9.0);
		double t159 = t41*y6*(4.0/3.0);
		double t160 = t45*y8*(4.0/3.0);
		double t161 = t46*y9*(4.0/3.0);
		double t162 = t156+t157+t158+t159+t160+t161;
		double t163 = t41*t91*t162*(1.0/6.0);
		double t164 = t89*t91*t104*(1.0/6.0);
		double t165 = t84*t91*t113*(1.0/6.0);
		double t193 = t92*y6*2.0;
		double t194 = t41*t104*t105*t113*(1.0/4.0);
		double t166 = t163+t164+t165-t193-t194;
		double t169 = t92*y0*2.0;
		double t170 = t75*t75;
		double t171 = t41*t105*t170*(1.0/4.0);
		double t172 = t41*y0*(4.0/3.0);
		double t173 = t45*y3*(4.0/3.0);
		double t174 = t46*y4*(4.0/3.0);
		double t175 = t79*t79;
		double t176 = t175*(2.0/9.0);
		double t177 = t124+t142+t172+t173+t174+t176;
		double t181 = t75*t79*t91*(1.0/3.0);
		double t182 = t41*t91*t177*(1.0/6.0);
		double t178 = t169+t171-t181-t182;
		double t183 = t92*y5*2.0;
		double t184 = t104*t104;
		double t185 = t41*t105*t184*(1.0/4.0);
		double t186 = t41*y5*(4.0/3.0);
		double t187 = t45*y1*(4.0/3.0);
		double t188 = t46*y8*(4.0/3.0);
		double t189 = t96*t96;
		double t190 = t189*(2.0/9.0);
		double t191 = t124+t144+t186+t187+t188+t190;
		double t196 = t84*t91*t104*(1.0/3.0);
		double t197 = t41*t91*t191*(1.0/6.0);
		double t192 = t183+t185-t196-t197;
		double t195 = t46*t72*t84*t117*(1.0/9.0);
		double t198 = t92*y7*2.0;
		double t199 = t113*t113;
		double t200 = t41*t105*t199*(1.0/4.0);
		double t201 = t41*y7*(4.0/3.0);
		double t202 = t46*y2*(4.0/3.0);
		double t203 = t45*y9*(4.0/3.0);
		double t204 = t109*t109;
		double t205 = t204*(2.0/9.0);
		double t206 = t142+t144+t201+t202+t203+t205;
		double t209 = t89*t91*t113*(1.0/3.0);
		double t210 = t41*t91*t206*(1.0/6.0);
		double t207 = t198+t200-t209-t210;
		double t208 = t45*t72*t89*t115*(1.0/9.0);
		double t211 = t41*t72*t100*t137*(1.0/9.0);
		return t74*(t74*t178+t137*(t41*t72*t79*(2.0/9.0)-t53*t75*t119*(1.0/9.0))-t41*t45*t72*t135*(1.0/9.0)+t45*t72*t79*t115*(1.0/9.0)+t41*t72*t84*t115*(1.0/9.0)-t41*t46*t72*t155*(1.0/9.0)+t46*t72*t79*t117*(1.0/9.0)+t41*t72*t89*t117*(1.0/9.0)-t41*t45*t75*t115*t119*(1.0/9.0)-t41*t46*t75*t117*t119*(1.0/9.0))+t122*(t122*t192+t115*(t45*t72*t96*(2.0/9.0)-t61*t104*t119*(1.0/9.0))-t41*t45*t72*t135*(1.0/9.0)-t45*t46*t72*t166*(1.0/9.0)+t46*t72*t96*t117*(1.0/9.0)+t45*t72*t100*t117*(1.0/9.0)+t45*t72*t84*t137*(1.0/9.0)+t41*t72*t96*t137*(1.0/9.0)-t45*t46*t104*t117*t119*(1.0/9.0)-t41*t45*t104*t119*t137*(1.0/9.0))+t140*(t140*t207+t117*(t46*t72*t109*(2.0/9.0)-t69*t113*t119*(1.0/9.0))-t41*t46*t72*t155*(1.0/9.0)-t45*t46*t72*t166*(1.0/9.0)+t46*t72*t100*t115*(1.0/9.0)+t45*t72*t109*t115*(1.0/9.0)+t46*t72*t89*t137*(1.0/9.0)+t41*t72*t109*t137*(1.0/9.0)-t45*t46*t113*t115*t119*(1.0/9.0)-t41*t46*t113*t119*t137*(1.0/9.0))+t41*t45*t72*(t195-t122*t135+t115*(t45*t72*t84*(2.0/9.0)-t61*t75*t119*(1.0/9.0))-t45*t46*t72*t155*(1.0/9.0)+t45*t72*t89*t117*(1.0/9.0)+t45*t72*t79*t137*(1.0/9.0)+t41*t72*t84*t137*(1.0/9.0)+t41*t45*t72*t178*(1.0/9.0)-t45*t46*t75*t117*t119*(1.0/9.0)-t41*t45*t75*t119*t137*(1.0/9.0))*(1.0/9.0)+t41*t45*t72*(t195-t74*t135+t137*(t41*t72*t84*(2.0/9.0)-t53*t104*t119*(1.0/9.0))+t45*t72*t84*t115*(1.0/9.0)+t41*t72*t96*t115*(1.0/9.0)-t41*t46*t72*t166*(1.0/9.0)+t41*t72*t100*t117*(1.0/9.0)+t41*t45*t72*t192*(1.0/9.0)-t41*t45*t104*t115*t119*(1.0/9.0)-t41*t46*t104*t117*t119*(1.0/9.0))*(1.0/9.0)+t41*t46*t72*(t208-t140*t155+t117*(t46*t72*t89*(2.0/9.0)-t69*t75*t119*(1.0/9.0))-t45*t46*t72*t135*(1.0/9.0)+t46*t72*t84*t115*(1.0/9.0)+t46*t72*t79*t137*(1.0/9.0)+t41*t46*t72*t178*(1.0/9.0)+t41*t72*t89*t137*(1.0/9.0)-t45*t46*t75*t115*t119*(1.0/9.0)-t41*t46*t75*t119*t137*(1.0/9.0))*(1.0/9.0)+t41*t46*t72*(t208-t74*t155+t137*(t41*t72*t89*(2.0/9.0)-t53*t113*t119*(1.0/9.0))-t41*t45*t72*t166*(1.0/9.0)+t46*t72*t89*t117*(1.0/9.0)+t41*t72*t100*t115*(1.0/9.0)+t41*t72*t109*t117*(1.0/9.0)+t41*t46*t72*t207*(1.0/9.0)-t41*t45*t113*t115*t119*(1.0/9.0)-t41*t46*t113*t117*t119*(1.0/9.0))*(1.0/9.0)+t45*t46*t72*(t211-t140*t166+t117*(t46*t72*t100*(2.0/9.0)-t69*t104*t119*(1.0/9.0))-t41*t46*t72*t135*(1.0/9.0)+t46*t72*t96*t115*(1.0/9.0)+t45*t72*t100*t115*(1.0/9.0)+t46*t72*t84*t137*(1.0/9.0)+t45*t46*t72*t192*(1.0/9.0)-t45*t46*t104*t115*t119*(1.0/9.0)-t41*t46*t104*t119*t137*(1.0/9.0))*(1.0/9.0)+t45*t46*t72*(t211-t122*t166+t115*(t45*t72*t100*(2.0/9.0)-t61*t113*t119*(1.0/9.0))-t41*t45*t72*t155*(1.0/9.0)+t46*t72*t100*t117*(1.0/9.0)+t45*t72*t89*t137*(1.0/9.0)+t45*t72*t109*t117*(1.0/9.0)+t45*t46*t72*t207*(1.0/9.0)-t45*t46*t113*t117*t119*(1.0/9.0)-t41*t45*t113*t119*t137*(1.0/9.0))*(1.0/9.0);
	}

	static inline double getLaplaceN1(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double y6 = y[6] ;
		double y7 = y[7] ;
		double y8 = y[8] ;
		double y9 = y[9] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t214 = x1*x1;
		double t215 = x0*x0;
		double t216 = x2*x2;
		double t219 = t215*y0*3.0;
		double t220 = t214*y5*3.0;
		double t221 = t216*y7*3.0;
		double t222 = x0*x1*y3*6.0;
		double t223 = x0*x2*y4*6.0;
		double t224 = x1*x2*y6*6.0;
		double t213 = t219+t220+t221+t222+t223+t224;
		double t227 = t214*y1*3.0;
		double t228 = t215*y3*3.0;
		double t229 = t216*y9*3.0;
		double t230 = x0*x1*y5*6.0;
		double t231 = x0*x2*y6*6.0;
		double t232 = x1*x2*y8*6.0;
		double t217 = t227+t228+t229+t230+t231+t232;
		double t235 = t215*y4*3.0;
		double t236 = t216*y2*3.0;
		double t237 = t214*y8*3.0;
		double t238 = x0*x1*y6*6.0;
		double t239 = x0*x2*y7*6.0;
		double t240 = x1*x2*y9*6.0;
		double t218 = t235+t236+t237+t238+t239+t240;
		double t225 = t213*t213;
		double t226 = t225*(1.0/9.0);
		double t233 = t217*t217;
		double t234 = t233*(1.0/9.0);
		double t241 = t218*t218;
		double t242 = t241*(1.0/9.0);
		double t243 = t226+t234+t242;
		double t244 = 1.0/t243;
		double t245 = t225*t244*(1.0/9.0);
		double t246 = t245-1.0;
		double t248 = x0*y3*6.0;
		double t249 = x1*y5*6.0;
		double t250 = x2*y6*6.0;
		double t251 = t248+t249+t250;
		double t252 = x0*y0*6.0;
		double t253 = x1*y3*6.0;
		double t254 = x2*y4*6.0;
		double t255 = t252+t253+t254;
		double t256 = t213*t255*(2.0/9.0);
		double t257 = t217*t251*(2.0/9.0);
		double t258 = x0*y4*6.0;
		double t259 = x1*y6*6.0;
		double t260 = x2*y7*6.0;
		double t261 = t258+t259+t260;
		double t262 = t218*t261*(2.0/9.0);
		double t247 = t256+t257+t262;
		double t263 = 1.0/pow(t243,3.0/2.0);
		double t264 = 1.0/sqrt(t243);
		double t265 = x1*y1*6.0;
		double t266 = x0*y5*6.0;
		double t267 = x2*y8*6.0;
		double t268 = t265+t266+t267;
		double t269 = x0*y6*6.0;
		double t270 = x1*y8*6.0;
		double t271 = x2*y9*6.0;
		double t272 = t269+t270+t271;
		double t273 = t213*t251*(2.0/9.0);
		double t274 = t217*t268*(2.0/9.0);
		double t275 = t218*t272*(2.0/9.0);
		double t276 = t273+t274+t275;
		double t277 = 1.0/pow(t243,5.0/2.0);
		double t278 = x2*y2*6.0;
		double t279 = x0*y7*6.0;
		double t280 = x1*y9*6.0;
		double t281 = t278+t279+t280;
		double t282 = t213*t261*(2.0/9.0);
		double t283 = t218*t281*(2.0/9.0);
		double t284 = t217*t272*(2.0/9.0);
		double t285 = t282+t283+t284;
		double t286 = t264*t268*(1.0/3.0);
		double t290 = t217*t263*t276*(1.0/6.0);
		double t287 = t286-t290;
		double t288 = t264*t272*(1.0/3.0);
		double t292 = t217*t263*t285*(1.0/6.0);
		double t289 = t288-t292;
		double t291 = 1.0/(t243*t243);
		double t293 = t233*t244*(1.0/9.0);
		double t294 = t293-1.0;
		double t295 = t251*t251;
		double t296 = t295*(2.0/9.0);
		double t297 = t251*t255*(2.0/9.0);
		double t298 = t251*t268*(2.0/9.0);
		double t299 = t261*t272*(2.0/9.0);
		double t300 = t213*y3*(4.0/3.0);
		double t301 = t217*y5*(4.0/3.0);
		double t302 = t218*y6*(4.0/3.0);
		double t303 = t297+t298+t299+t300+t301+t302;
		double t304 = t217*t263*t303*(1.0/6.0);
		double t305 = t247*t263*t268*(1.0/6.0);
		double t306 = t251*t263*t276*(1.0/6.0);
		double t339 = t264*y5*2.0;
		double t340 = t217*t247*t276*t277*(1.0/4.0);
		double t307 = t304+t305+t306-t339-t340;
		double t308 = t251*t264*(1.0/3.0);
		double t310 = t217*t247*t263*(1.0/6.0);
		double t309 = t308-t310;
		double t311 = t241*t244*(1.0/9.0);
		double t312 = t311-1.0;
		double t313 = t261*t261;
		double t314 = t313*(2.0/9.0);
		double t315 = t272*t272;
		double t316 = t315*(2.0/9.0);
		double t317 = t255*t261*(2.0/9.0);
		double t318 = t261*t281*(2.0/9.0);
		double t319 = t251*t272*(2.0/9.0);
		double t320 = t213*y4*(4.0/3.0);
		double t321 = t217*y6*(4.0/3.0);
		double t322 = t218*y7*(4.0/3.0);
		double t323 = t317+t318+t319+t320+t321+t322;
		double t324 = t217*t263*t323*(1.0/6.0);
		double t325 = t247*t263*t272*(1.0/6.0);
		double t326 = t251*t263*t285*(1.0/6.0);
		double t351 = t264*y6*2.0;
		double t352 = t217*t247*t277*t285*(1.0/4.0);
		double t327 = t324+t325+t326-t351-t352;
		double t328 = t251*t261*(2.0/9.0);
		double t329 = t268*t272*(2.0/9.0);
		double t330 = t272*t281*(2.0/9.0);
		double t331 = t213*y6*(4.0/3.0);
		double t332 = t217*y8*(4.0/3.0);
		double t333 = t218*y9*(4.0/3.0);
		double t334 = t328+t329+t330+t331+t332+t333;
		double t335 = t217*t263*t334*(1.0/6.0);
		double t336 = t263*t268*t285*(1.0/6.0);
		double t337 = t263*t272*t276*(1.0/6.0);
		double t365 = t264*y8*2.0;
		double t366 = t217*t276*t277*t285*(1.0/4.0);
		double t338 = t335+t336+t337-t365-t366;
		double t341 = t264*y3*2.0;
		double t342 = t247*t247;
		double t343 = t217*t277*t342*(1.0/4.0);
		double t344 = t213*y0*(4.0/3.0);
		double t345 = t217*y3*(4.0/3.0);
		double t346 = t218*y4*(4.0/3.0);
		double t347 = t255*t255;
		double t348 = t347*(2.0/9.0);
		double t349 = t296+t314+t344+t345+t346+t348;
		double t353 = t247*t251*t263*(1.0/3.0);
		double t354 = t217*t263*t349*(1.0/6.0);
		double t350 = t341+t343-t353-t354;
		double t355 = t264*y1*2.0;
		double t356 = t276*t276;
		double t357 = t217*t277*t356*(1.0/4.0);
		double t358 = t213*y5*(4.0/3.0);
		double t359 = t217*y1*(4.0/3.0);
		double t360 = t218*y8*(4.0/3.0);
		double t361 = t268*t268;
		double t362 = t361*(2.0/9.0);
		double t363 = t296+t316+t358+t359+t360+t362;
		double t368 = t263*t268*t276*(1.0/3.0);
		double t369 = t217*t263*t363*(1.0/6.0);
		double t364 = t355+t357-t368-t369;
		double t367 = t218*t244*t251*t289*(1.0/9.0);
		double t370 = t264*y9*2.0;
		double t371 = t285*t285;
		double t372 = t217*t277*t371*(1.0/4.0);
		double t373 = t213*y7*(4.0/3.0);
		double t374 = t218*y2*(4.0/3.0);
		double t375 = t217*y9*(4.0/3.0);
		double t376 = t281*t281;
		double t377 = t376*(2.0/9.0);
		double t378 = t314+t316+t373+t374+t375+t377;
		double t381 = t263*t272*t285*(1.0/3.0);
		double t382 = t217*t263*t378*(1.0/6.0);
		double t379 = t370+t372-t381-t382;
		double t380 = t217*t244*t261*t287*(1.0/9.0);
		double t383 = t213*t244*t272*t309*(1.0/9.0);
		return t246*(t246*t350+t309*(t213*t244*t255*(2.0/9.0)-t225*t247*t291*(1.0/9.0))-t213*t217*t244*t307*(1.0/9.0)+t213*t244*t251*t287*(1.0/9.0)-t213*t218*t244*t327*(1.0/9.0)+t217*t244*t255*t287*(1.0/9.0)+t218*t244*t255*t289*(1.0/9.0)+t213*t244*t261*t289*(1.0/9.0)-t213*t217*t247*t287*t291*(1.0/9.0)-t213*t218*t247*t289*t291*(1.0/9.0))+t294*(t294*t364+t287*(t217*t244*t268*(2.0/9.0)-t233*t276*t291*(1.0/9.0))-t213*t217*t244*t307*(1.0/9.0)-t217*t218*t244*t338*(1.0/9.0)+t218*t244*t268*t289*(1.0/9.0)+t217*t244*t251*t309*(1.0/9.0)+t217*t244*t272*t289*(1.0/9.0)+t213*t244*t268*t309*(1.0/9.0)-t217*t218*t276*t289*t291*(1.0/9.0)-t213*t217*t276*t291*t309*(1.0/9.0))+t312*(t312*t379+t289*(t218*t244*t281*(2.0/9.0)-t241*t285*t291*(1.0/9.0))-t213*t218*t244*t327*(1.0/9.0)-t217*t218*t244*t338*(1.0/9.0)+t218*t244*t272*t287*(1.0/9.0)+t217*t244*t281*t287*(1.0/9.0)+t218*t244*t261*t309*(1.0/9.0)+t213*t244*t281*t309*(1.0/9.0)-t217*t218*t285*t287*t291*(1.0/9.0)-t213*t218*t285*t291*t309*(1.0/9.0))+t213*t217*t244*(t367-t294*t307+t287*(t217*t244*t251*(2.0/9.0)-t233*t247*t291*(1.0/9.0))-t217*t218*t244*t327*(1.0/9.0)+t217*t244*t261*t289*(1.0/9.0)+t213*t244*t251*t309*(1.0/9.0)+t213*t217*t244*t350*(1.0/9.0)+t217*t244*t255*t309*(1.0/9.0)-t217*t218*t247*t289*t291*(1.0/9.0)-t213*t217*t247*t291*t309*(1.0/9.0))*(1.0/9.0)+t213*t217*t244*(t367-t246*t307+t309*(t213*t244*t251*(2.0/9.0)-t225*t276*t291*(1.0/9.0))+t217*t244*t251*t287*(1.0/9.0)+t213*t244*t268*t287*(1.0/9.0)-t213*t218*t244*t338*(1.0/9.0)+t213*t244*t272*t289*(1.0/9.0)+t213*t217*t244*t364*(1.0/9.0)-t213*t217*t276*t287*t291*(1.0/9.0)-t213*t218*t276*t289*t291*(1.0/9.0))*(1.0/9.0)+t213*t218*t244*(t380-t312*t327+t289*(t218*t244*t261*(2.0/9.0)-t241*t247*t291*(1.0/9.0))-t217*t218*t244*t307*(1.0/9.0)+t218*t244*t251*t287*(1.0/9.0)+t213*t218*t244*t350*(1.0/9.0)+t218*t244*t255*t309*(1.0/9.0)+t213*t244*t261*t309*(1.0/9.0)-t217*t218*t247*t287*t291*(1.0/9.0)-t213*t218*t247*t291*t309*(1.0/9.0))*(1.0/9.0)+t213*t218*t244*(t380-t246*t327+t309*(t213*t244*t261*(2.0/9.0)-t225*t285*t291*(1.0/9.0))-t213*t217*t244*t338*(1.0/9.0)+t218*t244*t261*t289*(1.0/9.0)+t213*t244*t272*t287*(1.0/9.0)+t213*t244*t281*t289*(1.0/9.0)+t213*t218*t244*t379*(1.0/9.0)-t213*t217*t285*t287*t291*(1.0/9.0)-t213*t218*t285*t289*t291*(1.0/9.0))*(1.0/9.0)+t217*t218*t244*(t383-t312*t338+t289*(t218*t244*t272*(2.0/9.0)-t241*t276*t291*(1.0/9.0))-t213*t218*t244*t307*(1.0/9.0)+t218*t244*t268*t287*(1.0/9.0)+t217*t244*t272*t287*(1.0/9.0)+t218*t244*t251*t309*(1.0/9.0)+t217*t218*t244*t364*(1.0/9.0)-t217*t218*t276*t287*t291*(1.0/9.0)-t213*t218*t276*t291*t309*(1.0/9.0))*(1.0/9.0)+t217*t218*t244*(t383-t294*t338+t287*(t217*t244*t272*(2.0/9.0)-t233*t285*t291*(1.0/9.0))-t213*t217*t244*t327*(1.0/9.0)+t218*t244*t272*t289*(1.0/9.0)+t217*t244*t261*t309*(1.0/9.0)+t217*t244*t281*t289*(1.0/9.0)+t217*t218*t244*t379*(1.0/9.0)-t217*t218*t285*t289*t291*(1.0/9.0)-t213*t217*t285*t291*t309*(1.0/9.0))*(1.0/9.0);
	}

	static inline double getLaplaceN2(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double y6 = y[6] ;
		double y7 = y[7] ;
		double y8 = y[8] ;
		double y9 = y[9] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t386 = x1*x1;
		double t387 = x0*x0;
		double t388 = x2*x2;
		double t391 = t387*y0*3.0;
		double t392 = t386*y5*3.0;
		double t393 = t388*y7*3.0;
		double t394 = x0*x1*y3*6.0;
		double t395 = x0*x2*y4*6.0;
		double t396 = x1*x2*y6*6.0;
		double t385 = t391+t392+t393+t394+t395+t396;
		double t399 = t386*y1*3.0;
		double t400 = t387*y3*3.0;
		double t401 = t388*y9*3.0;
		double t402 = x0*x1*y5*6.0;
		double t403 = x0*x2*y6*6.0;
		double t404 = x1*x2*y8*6.0;
		double t389 = t399+t400+t401+t402+t403+t404;
		double t407 = t387*y4*3.0;
		double t408 = t388*y2*3.0;
		double t409 = t386*y8*3.0;
		double t410 = x0*x1*y6*6.0;
		double t411 = x0*x2*y7*6.0;
		double t412 = x1*x2*y9*6.0;
		double t390 = t407+t408+t409+t410+t411+t412;
		double t397 = t385*t385;
		double t398 = t397*(1.0/9.0);
		double t405 = t389*t389;
		double t406 = t405*(1.0/9.0);
		double t413 = t390*t390;
		double t414 = t413*(1.0/9.0);
		double t415 = t398+t406+t414;
		double t416 = 1.0/t415;
		double t417 = t397*t416*(1.0/9.0);
		double t418 = t417-1.0;
		double t420 = x0*y4*6.0;
		double t421 = x1*y6*6.0;
		double t422 = x2*y7*6.0;
		double t423 = t420+t421+t422;
		double t424 = x0*y0*6.0;
		double t425 = x1*y3*6.0;
		double t426 = x2*y4*6.0;
		double t427 = t424+t425+t426;
		double t428 = t385*t427*(2.0/9.0);
		double t429 = x0*y3*6.0;
		double t430 = x1*y5*6.0;
		double t431 = x2*y6*6.0;
		double t432 = t429+t430+t431;
		double t433 = t389*t432*(2.0/9.0);
		double t434 = t390*t423*(2.0/9.0);
		double t419 = t428+t433+t434;
		double t435 = 1.0/pow(t415,3.0/2.0);
		double t436 = 1.0/sqrt(t415);
		double t437 = x0*y6*6.0;
		double t438 = x1*y8*6.0;
		double t439 = x2*y9*6.0;
		double t440 = t437+t438+t439;
		double t441 = x1*y1*6.0;
		double t442 = x0*y5*6.0;
		double t443 = x2*y8*6.0;
		double t444 = t441+t442+t443;
		double t445 = t385*t432*(2.0/9.0);
		double t446 = t389*t444*(2.0/9.0);
		double t447 = t390*t440*(2.0/9.0);
		double t448 = t445+t446+t447;
		double t449 = 1.0/pow(t415,5.0/2.0);
		double t450 = x2*y2*6.0;
		double t451 = x0*y7*6.0;
		double t452 = x1*y9*6.0;
		double t453 = t450+t451+t452;
		double t454 = t385*t423*(2.0/9.0);
		double t455 = t390*t453*(2.0/9.0);
		double t456 = t389*t440*(2.0/9.0);
		double t457 = t454+t455+t456;
		double t458 = t436*t440*(1.0/3.0);
		double t462 = t390*t435*t448*(1.0/6.0);
		double t459 = t458-t462;
		double t460 = t436*t453*(1.0/3.0);
		double t464 = t390*t435*t457*(1.0/6.0);
		double t461 = t460-t464;
		double t463 = 1.0/(t415*t415);
		double t465 = t405*t416*(1.0/9.0);
		double t466 = t465-1.0;
		double t467 = t432*t432;
		double t468 = t467*(2.0/9.0);
		double t469 = t427*t432*(2.0/9.0);
		double t470 = t432*t444*(2.0/9.0);
		double t471 = t423*t440*(2.0/9.0);
		double t472 = t385*y3*(4.0/3.0);
		double t473 = t389*y5*(4.0/3.0);
		double t474 = t390*y6*(4.0/3.0);
		double t475 = t469+t470+t471+t472+t473+t474;
		double t476 = t390*t435*t475*(1.0/6.0);
		double t477 = t419*t435*t440*(1.0/6.0);
		double t478 = t423*t435*t448*(1.0/6.0);
		double t511 = t436*y6*2.0;
		double t512 = t390*t419*t448*t449*(1.0/4.0);
		double t479 = t476+t477+t478-t511-t512;
		double t480 = t423*t436*(1.0/3.0);
		double t482 = t390*t419*t435*(1.0/6.0);
		double t481 = t480-t482;
		double t483 = t413*t416*(1.0/9.0);
		double t484 = t483-1.0;
		double t485 = t423*t423;
		double t486 = t485*(2.0/9.0);
		double t487 = t440*t440;
		double t488 = t487*(2.0/9.0);
		double t489 = t423*t427*(2.0/9.0);
		double t490 = t423*t453*(2.0/9.0);
		double t491 = t432*t440*(2.0/9.0);
		double t492 = t385*y4*(4.0/3.0);
		double t493 = t389*y6*(4.0/3.0);
		double t494 = t390*y7*(4.0/3.0);
		double t495 = t489+t490+t491+t492+t493+t494;
		double t496 = t390*t435*t495*(1.0/6.0);
		double t497 = t419*t435*t453*(1.0/6.0);
		double t498 = t423*t435*t457*(1.0/6.0);
		double t523 = t436*y7*2.0;
		double t524 = t390*t419*t449*t457*(1.0/4.0);
		double t499 = t496+t497+t498-t523-t524;
		double t500 = t423*t432*(2.0/9.0);
		double t501 = t440*t444*(2.0/9.0);
		double t502 = t440*t453*(2.0/9.0);
		double t503 = t385*y6*(4.0/3.0);
		double t504 = t389*y8*(4.0/3.0);
		double t505 = t390*y9*(4.0/3.0);
		double t506 = t500+t501+t502+t503+t504+t505;
		double t507 = t390*t435*t506*(1.0/6.0);
		double t508 = t435*t448*t453*(1.0/6.0);
		double t509 = t435*t440*t457*(1.0/6.0);
		double t537 = t436*y9*2.0;
		double t538 = t390*t448*t449*t457*(1.0/4.0);
		double t510 = t507+t508+t509-t537-t538;
		double t513 = t436*y4*2.0;
		double t514 = t419*t419;
		double t515 = t390*t449*t514*(1.0/4.0);
		double t516 = t385*y0*(4.0/3.0);
		double t517 = t389*y3*(4.0/3.0);
		double t518 = t390*y4*(4.0/3.0);
		double t519 = t427*t427;
		double t520 = t519*(2.0/9.0);
		double t521 = t468+t486+t516+t517+t518+t520;
		double t525 = t419*t423*t435*(1.0/3.0);
		double t526 = t390*t435*t521*(1.0/6.0);
		double t522 = t513+t515-t525-t526;
		double t527 = t436*y8*2.0;
		double t528 = t448*t448;
		double t529 = t390*t449*t528*(1.0/4.0);
		double t530 = t385*y5*(4.0/3.0);
		double t531 = t389*y1*(4.0/3.0);
		double t532 = t390*y8*(4.0/3.0);
		double t533 = t444*t444;
		double t534 = t533*(2.0/9.0);
		double t535 = t468+t488+t530+t531+t532+t534;
		double t540 = t435*t440*t448*(1.0/3.0);
		double t541 = t390*t435*t535*(1.0/6.0);
		double t536 = t527+t529-t540-t541;
		double t539 = t390*t416*t432*t461*(1.0/9.0);
		double t542 = t436*y2*2.0;
		double t543 = t457*t457;
		double t544 = t390*t449*t543*(1.0/4.0);
		double t545 = t385*y7*(4.0/3.0);
		double t546 = t390*y2*(4.0/3.0);
		double t547 = t389*y9*(4.0/3.0);
		double t548 = t453*t453;
		double t549 = t548*(2.0/9.0);
		double t550 = t486+t488+t545+t546+t547+t549;
		double t553 = t435*t453*t457*(1.0/3.0);
		double t554 = t390*t435*t550*(1.0/6.0);
		double t551 = t542+t544-t553-t554;
		double t552 = t389*t416*t423*t459*(1.0/9.0);
		double t555 = t385*t416*t440*t481*(1.0/9.0);
		return t418*(t418*t522+t481*(t385*t416*t427*(2.0/9.0)-t397*t419*t463*(1.0/9.0))-t385*t389*t416*t479*(1.0/9.0)+t385*t416*t423*t461*(1.0/9.0)-t385*t390*t416*t499*(1.0/9.0)+t389*t416*t427*t459*(1.0/9.0)+t385*t416*t432*t459*(1.0/9.0)+t390*t416*t427*t461*(1.0/9.0)-t385*t389*t419*t459*t463*(1.0/9.0)-t385*t390*t419*t461*t463*(1.0/9.0))+t466*(t466*t536+t459*(t389*t416*t444*(2.0/9.0)-t405*t448*t463*(1.0/9.0))-t385*t389*t416*t479*(1.0/9.0)-t389*t390*t416*t510*(1.0/9.0)+t389*t416*t440*t461*(1.0/9.0)+t390*t416*t444*t461*(1.0/9.0)+t389*t416*t432*t481*(1.0/9.0)+t385*t416*t444*t481*(1.0/9.0)-t389*t390*t448*t461*t463*(1.0/9.0)-t385*t389*t448*t463*t481*(1.0/9.0))+t484*(t484*t551+t461*(t390*t416*t453*(2.0/9.0)-t413*t457*t463*(1.0/9.0))-t385*t390*t416*t499*(1.0/9.0)-t389*t390*t416*t510*(1.0/9.0)+t390*t416*t440*t459*(1.0/9.0)+t390*t416*t423*t481*(1.0/9.0)+t389*t416*t453*t459*(1.0/9.0)+t385*t416*t453*t481*(1.0/9.0)-t389*t390*t457*t459*t463*(1.0/9.0)-t385*t390*t457*t463*t481*(1.0/9.0))+t385*t389*t416*(t539-t466*t479+t459*(t389*t416*t432*(2.0/9.0)-t405*t419*t463*(1.0/9.0))+t389*t416*t423*t461*(1.0/9.0)-t389*t390*t416*t499*(1.0/9.0)+t385*t389*t416*t522*(1.0/9.0)+t389*t416*t427*t481*(1.0/9.0)+t385*t416*t432*t481*(1.0/9.0)-t389*t390*t419*t461*t463*(1.0/9.0)-t385*t389*t419*t463*t481*(1.0/9.0))*(1.0/9.0)+t385*t389*t416*(t539-t418*t479+t481*(t385*t416*t432*(2.0/9.0)-t397*t448*t463*(1.0/9.0))+t389*t416*t432*t459*(1.0/9.0)-t385*t390*t416*t510*(1.0/9.0)+t385*t416*t440*t461*(1.0/9.0)+t385*t416*t444*t459*(1.0/9.0)+t385*t389*t416*t536*(1.0/9.0)-t385*t389*t448*t459*t463*(1.0/9.0)-t385*t390*t448*t461*t463*(1.0/9.0))*(1.0/9.0)+t385*t390*t416*(t552-t484*t499+t461*(t390*t416*t423*(2.0/9.0)-t413*t419*t463*(1.0/9.0))-t389*t390*t416*t479*(1.0/9.0)+t390*t416*t432*t459*(1.0/9.0)+t385*t416*t423*t481*(1.0/9.0)+t385*t390*t416*t522*(1.0/9.0)+t390*t416*t427*t481*(1.0/9.0)-t389*t390*t419*t459*t463*(1.0/9.0)-t385*t390*t419*t463*t481*(1.0/9.0))*(1.0/9.0)+t385*t390*t416*(t552-t418*t499+t481*(t385*t416*t423*(2.0/9.0)-t397*t457*t463*(1.0/9.0))+t390*t416*t423*t461*(1.0/9.0)-t385*t389*t416*t510*(1.0/9.0)+t385*t416*t440*t459*(1.0/9.0)+t385*t416*t453*t461*(1.0/9.0)+t385*t390*t416*t551*(1.0/9.0)-t385*t389*t457*t459*t463*(1.0/9.0)-t385*t390*t457*t461*t463*(1.0/9.0))*(1.0/9.0)+t389*t390*t416*(t555-t484*t510+t461*(t390*t416*t440*(2.0/9.0)-t413*t448*t463*(1.0/9.0))-t385*t390*t416*t479*(1.0/9.0)+t389*t416*t440*t459*(1.0/9.0)+t390*t416*t444*t459*(1.0/9.0)+t390*t416*t432*t481*(1.0/9.0)+t389*t390*t416*t536*(1.0/9.0)-t389*t390*t448*t459*t463*(1.0/9.0)-t385*t390*t448*t463*t481*(1.0/9.0))*(1.0/9.0)+t389*t390*t416*(t555-t466*t510+t459*(t389*t416*t440*(2.0/9.0)-t405*t457*t463*(1.0/9.0))-t385*t389*t416*t499*(1.0/9.0)+t390*t416*t440*t461*(1.0/9.0)+t389*t416*t423*t481*(1.0/9.0)+t389*t416*t453*t461*(1.0/9.0)+t389*t390*t416*t551*(1.0/9.0)-t389*t390*t457*t461*t463*(1.0/9.0)-t385*t389*t457*t463*t481*(1.0/9.0))*(1.0/9.0);
	}

	static inline double getGaussianCurvature(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double y6 = y[6] ;
		double y7 = y[7] ;
		double y8 = y[8] ;
		double y9 = y[9] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t3 = x1*x1;
		double t4 = x0*x0;
		double t5 = x2*x2;
		double t8 = t4*y0*3.0;
		double t9 = t3*y5*3.0;
		double t10 = t5*y7*3.0;
		double t11 = x0*x1*y3*6.0;
		double t12 = x0*x2*y4*6.0;
		double t13 = x1*x2*y6*6.0;
		double t2 = t8+t9+t10+t11+t12+t13;
		double t16 = t3*y1*3.0;
		double t17 = t4*y3*3.0;
		double t18 = t5*y9*3.0;
		double t19 = x0*x1*y5*6.0;
		double t20 = x0*x2*y6*6.0;
		double t21 = x1*x2*y8*6.0;
		double t6 = t16+t17+t18+t19+t20+t21;
		double t24 = t4*y4*3.0;
		double t25 = t5*y2*3.0;
		double t26 = t3*y8*3.0;
		double t27 = x0*x1*y6*6.0;
		double t28 = x0*x2*y7*6.0;
		double t29 = x1*x2*y9*6.0;
		double t7 = t24+t25+t26+t27+t28+t29;
		double t14 = t2*t2;
		double t15 = t14*(1.0/9.0);
		double t22 = t6*t6;
		double t23 = t22*(1.0/9.0);
		double t30 = t7*t7;
		double t31 = t30*(1.0/9.0);
		double t32 = t15+t23+t31;
		double t33 = x0*y0*6.0;
		double t34 = x1*y3*6.0;
		double t35 = x2*y4*6.0;
		double t36 = t33+t34+t35;
		double t37 = 1.0/t32;
		double t38 = 1.0/sqrt(t32);
		double t39 = x0*y3*6.0;
		double t40 = x1*y5*6.0;
		double t41 = x2*y6*6.0;
		double t42 = t39+t40+t41;
		double t43 = x1*y1*6.0;
		double t44 = x0*y5*6.0;
		double t45 = x2*y8*6.0;
		double t46 = t43+t44+t45;
		double t47 = 1.0/pow(t32,3.0/2.0);
		double t48 = x0*y4*6.0;
		double t49 = x1*y6*6.0;
		double t50 = x2*y7*6.0;
		double t51 = t48+t49+t50;
		double t52 = x2*y2*6.0;
		double t53 = x0*y7*6.0;
		double t54 = x1*y9*6.0;
		double t55 = t52+t53+t54;
		double t56 = x0*y6*6.0;
		double t57 = x1*y8*6.0;
		double t58 = x2*y9*6.0;
		double t59 = t56+t57+t58;
		double t60 = t2*t36*(2.0/9.0);
		double t61 = t6*t42*(2.0/9.0);
		double t62 = t7*t51*(2.0/9.0);
		double t63 = t60+t61+t62;
		double t64 = t38*t42*(1.0/3.0);
		double t65 = t2*t42*(2.0/9.0);
		double t66 = t6*t46*(2.0/9.0);
		double t67 = t7*t59*(2.0/9.0);
		double t68 = t65+t66+t67;
		double t69 = t38*t51*(1.0/3.0);
		double t70 = t2*t51*(2.0/9.0);
		double t71 = t7*t55*(2.0/9.0);
		double t72 = t6*t59*(2.0/9.0);
		double t73 = t70+t71+t72;
		double t74 = t38*t59*(1.0/3.0);
		double t76 = t14*t37*(1.0/9.0);
		double t77 = t76-1.0;
		double t78 = t36*t38*(1.0/3.0);
		double t79 = t2*t47*t63*(1.0/6.0);
		double t80 = t78-t79;
		double t81 = t77*t80;
		double t82 = t2*t47*t68*(1.0/6.0);
		double t83 = t64-t82;
		double t84 = t2*t6*t37*t83*(1.0/9.0);
		double t85 = t2*t47*t73*(1.0/6.0);
		double t86 = t69-t85;
		double t87 = t2*t7*t37*t86*(1.0/9.0);
		double t89 = t22*t37*(1.0/9.0);
		double t90 = t89-1.0;
		double t92 = t30*t37*(1.0/9.0);
		double t93 = t92-1.0;
		double t95 = t6*t47*t63*(1.0/6.0);
		double t96 = t64-t95;
		double t97 = t38*t46*(1.0/3.0);
		double t98 = t6*t47*t68*(1.0/6.0);
		double t99 = t97-t98;
		double t100 = t6*t47*t73*(1.0/6.0);
		double t101 = t74-t100;
		double t103 = t90*t99;
		double t104 = t2*t6*t37*t96*(1.0/9.0);
		double t105 = t6*t7*t37*t101*(1.0/9.0);
		double t107 = t7*t47*t63*(1.0/6.0);
		double t108 = t69-t107;
		double t109 = t7*t47*t68*(1.0/6.0);
		double t110 = t74-t109;
		double t111 = t38*t55*(1.0/3.0);
		double t112 = t7*t47*t73*(1.0/6.0);
		double t113 = t111-t112;
		double t117 = t93*t113;
		double t118 = t2*t7*t37*t108*(1.0/9.0);
		double t119 = t6*t7*t37*t110*(1.0/9.0);
		double t75 = t81+t84+t87+t103+t104+t105+t117+t118+t119;
		double t88 = t81+t84+t87;
		double t91 = t83*t90+t2*t6*t37*t80*(1.0/9.0)+t6*t7*t37*t86*(1.0/9.0);
		double t94 = t86*t93+t2*t7*t37*t80*(1.0/9.0)+t6*t7*t37*t83*(1.0/9.0);
		double t102 = t77*t96+t2*t6*t37*t99*(1.0/9.0)+t2*t7*t37*t101*(1.0/9.0);
		double t106 = t103+t104+t105;
		double t114 = t77*t108+t2*t6*t37*t110*(1.0/9.0)+t2*t7*t37*t113*(1.0/9.0);
		double t115 = t93*t101+t2*t7*t37*t96*(1.0/9.0)+t6*t7*t37*t99*(1.0/9.0);
		double t116 = t90*t110+t2*t6*t37*t108*(1.0/9.0)+t6*t7*t37*t113*(1.0/9.0);
		double t120 = t117+t118+t119;
		return (t75*t75)*(1.0/2.0)-(t88*t88)*(1.0/2.0)-(t91*t91)*(1.0/2.0)-(t94*t94)*(1.0/2.0)-(t102*t102)*(1.0/2.0)-(t106*t106)*(1.0/2.0)-(t114*t114)*(1.0/2.0)-(t115*t115)*(1.0/2.0)-(t116*t116)*(1.0/2.0)-(t120*t120)*(1.0/2.0);
	}

	static inline double getMeanCurvature(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double y6 = y[6] ;
		double y7 = y[7] ;
		double y8 = y[8] ;
		double y9 = y[9] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t1017 = x1*x1;
		double t1018 = x0*x0;
		double t1019 = x2*x2;
		double t1022 = t1018*y0*3.0;
		double t1023 = t1017*y5*3.0;
		double t1024 = t1019*y7*3.0;
		double t1025 = x0*x1*y3*6.0;
		double t1026 = x0*x2*y4*6.0;
		double t1027 = x1*x2*y6*6.0;
		double t1016 = t1022+t1023+t1024+t1025+t1026+t1027;
		double t1030 = t1017*y1*3.0;
		double t1031 = t1018*y3*3.0;
		double t1032 = t1019*y9*3.0;
		double t1033 = x0*x1*y5*6.0;
		double t1034 = x0*x2*y6*6.0;
		double t1035 = x1*x2*y8*6.0;
		double t1020 = t1030+t1031+t1032+t1033+t1034+t1035;
		double t1038 = t1018*y4*3.0;
		double t1039 = t1019*y2*3.0;
		double t1040 = t1017*y8*3.0;
		double t1041 = x0*x1*y6*6.0;
		double t1042 = x0*x2*y7*6.0;
		double t1043 = x1*x2*y9*6.0;
		double t1021 = t1038+t1039+t1040+t1041+t1042+t1043;
		double t1028 = t1016*t1016;
		double t1029 = t1028*(1.0/9.0);
		double t1036 = t1020*t1020;
		double t1037 = t1036*(1.0/9.0);
		double t1044 = t1021*t1021;
		double t1045 = t1044*(1.0/9.0);
		double t1046 = t1029+t1037+t1045;
		double t1047 = x0*y0*6.0;
		double t1048 = x1*y3*6.0;
		double t1049 = x2*y4*6.0;
		double t1050 = t1047+t1048+t1049;
		double t1051 = 1.0/t1046;
		double t1052 = 1.0/sqrt(t1046);
		double t1053 = x0*y3*6.0;
		double t1054 = x1*y5*6.0;
		double t1055 = x2*y6*6.0;
		double t1056 = t1053+t1054+t1055;
		double t1057 = x1*y1*6.0;
		double t1058 = x0*y5*6.0;
		double t1059 = x2*y8*6.0;
		double t1060 = t1057+t1058+t1059;
		double t1061 = 1.0/pow(t1046,3.0/2.0);
		double t1062 = x0*y4*6.0;
		double t1063 = x1*y6*6.0;
		double t1064 = x2*y7*6.0;
		double t1065 = t1062+t1063+t1064;
		double t1066 = x2*y2*6.0;
		double t1067 = x0*y7*6.0;
		double t1068 = x1*y9*6.0;
		double t1069 = t1066+t1067+t1068;
		double t1070 = x0*y6*6.0;
		double t1071 = x1*y8*6.0;
		double t1072 = x2*y9*6.0;
		double t1073 = t1070+t1071+t1072;
		double t1074 = t1016*t1050*(2.0/9.0);
		double t1075 = t1020*t1056*(2.0/9.0);
		double t1076 = t1021*t1065*(2.0/9.0);
		double t1077 = t1074+t1075+t1076;
		double t1078 = t1052*t1056*(1.0/3.0);
		double t1079 = t1016*t1056*(2.0/9.0);
		double t1080 = t1020*t1060*(2.0/9.0);
		double t1081 = t1021*t1073*(2.0/9.0);
		double t1082 = t1079+t1080+t1081;
		double t1083 = t1052*t1065*(1.0/3.0);
		double t1084 = t1016*t1065*(2.0/9.0);
		double t1085 = t1021*t1069*(2.0/9.0);
		double t1086 = t1020*t1073*(2.0/9.0);
		double t1087 = t1084+t1085+t1086;
		double t1088 = t1052*t1073*(1.0/3.0);
		return -(t1050*t1052*(1.0/3.0)-t1016*t1061*t1077*(1.0/6.0))*(t1028*t1051*(1.0/9.0)-1.0)-(t1052*t1060*(1.0/3.0)-t1020*t1061*t1082*(1.0/6.0))*(t1036*t1051*(1.0/9.0)-1.0)-(t1052*t1069*(1.0/3.0)-t1021*t1061*t1087*(1.0/6.0))*(t1044*t1051*(1.0/9.0)-1.0)-t1016*t1020*t1051*(t1078-t1020*t1061*t1077*(1.0/6.0))*(1.0/9.0)-t1016*t1020*t1051*(t1078-t1016*t1061*t1082*(1.0/6.0))*(1.0/9.0)-t1016*t1021*t1051*(t1083-t1021*t1061*t1077*(1.0/6.0))*(1.0/9.0)-t1016*t1021*t1051*(t1083-t1016*t1061*t1087*(1.0/6.0))*(1.0/9.0)-t1020*t1021*t1051*(t1088-t1021*t1061*t1082*(1.0/6.0))*(1.0/9.0)-t1020*t1021*t1051*(t1088-t1020*t1061*t1087*(1.0/6.0))*(1.0/9.0);
	}

private:
};



#endif /* QUBICGEOMETRYQUANTITIES_H_ */
