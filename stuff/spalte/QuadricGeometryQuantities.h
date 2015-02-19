/*
 * QuadricGeometryQuantities.h
 *
 *  Created on: 23.05.2014
 *      Author: sebastian
 */

#ifndef QUADRICGEOMETRYQUANTITIES_H_
#define QUADRICGEOMETRYQUANTITIES_H_

class QuadricGeometryQuantities
{
public:
	QuadricGeometryQuantities() {}

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
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t2 = x0*y0*2.0;
		double t3 = x1*y3*2.0;
		double t4 = x2*y4*2.0;
		double t5 = t2+t3+t4;
		double t6 = x1*y1*2.0+x0*y3*2.0+x2*y5*2.0;
		double t7 = x0*y4*2.0+x2*y2*2.0+x1*y5*2.0;
		return t5*1.0/sqrt((t5*t5)*(1.0/4.0)+(t6*t6)*(1.0/4.0)+(t7*t7)*(1.0/4.0))*(1.0/2.0);
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
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t9 = x0*y0*2.0+x1*y3*2.0+x2*y4*2.0;
		double t10 = x1*y1*2.0;
		double t11 = x0*y3*2.0;
		double t12 = x2*y5*2.0;
		double t13 = t10+t11+t12;
		double t14 = x0*y4*2.0+x2*y2*2.0+x1*y5*2.0;
		return t13*1.0/sqrt((t9*t9)*(1.0/4.0)+(t13*t13)*(1.0/4.0)+(t14*t14)*(1.0/4.0))*(1.0/2.0);
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
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t16 = x0*y0*2.0+x1*y3*2.0+x2*y4*2.0;
		double t17 = x1*y1*2.0+x0*y3*2.0+x2*y5*2.0;
		double t18 = x0*y4*2.0;
		double t19 = x2*y2*2.0;
		double t20 = x1*y5*2.0;
		double t21 = t18+t19+t20;
		return t21*1.0/sqrt((t16*t16)*(1.0/4.0)+(t17*t17)*(1.0/4.0)+(t21*t21)*(1.0/4.0))*(1.0/2.0);
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
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t24 = x0*y0*2.0;
		double t25 = x1*y3*2.0;
		double t26 = x2*y4*2.0;
		double t23 = t24+t25+t26;
		double t27 = t23*t23;
		double t31 = x1*y1*2.0;
		double t32 = x0*y3*2.0;
		double t33 = x2*y5*2.0;
		double t28 = t31+t32+t33;
		double t36 = x0*y4*2.0;
		double t37 = x2*y2*2.0;
		double t38 = x1*y5*2.0;
		double t29 = t36+t37+t38;
		double t30 = t27*(1.0/4.0);
		double t34 = t28*t28;
		double t35 = t34*(1.0/4.0);
		double t39 = t29*t29;
		double t40 = t39*(1.0/4.0);
		double t41 = t30+t35+t40;
		double t42 = 1.0/t41;
		double t43 = t27*t42*(1.0/4.0);
		double t44 = t43-1.0;
		double t45 = 1.0/pow(t41,3.0/2.0);
		double t46 = t23*y0;
		double t47 = t28*y3;
		double t48 = t29*y4;
		double t49 = t46+t47+t48;
		double t50 = 1.0/sqrt(t41);
		double t51 = t50*y3;
		double t52 = t28*y1;
		double t53 = t23*y3;
		double t54 = t29*y5;
		double t55 = t52+t53+t54;
		double t64 = t23*t45*t55*(1.0/4.0);
		double t56 = t51-t64;
		double t57 = t50*y4;
		double t58 = t23*y4;
		double t59 = t29*y2;
		double t60 = t28*y5;
		double t61 = t58+t59+t60;
		double t66 = t23*t45*t61*(1.0/4.0);
		double t62 = t57-t66;
		double t63 = 1.0/pow(t41,5.0/2.0);
		double t65 = 1.0/(t41*t41);
		double t67 = t34*t42*(1.0/4.0);
		double t68 = t67-1.0;
		double t69 = y3*y3;
		double t70 = t69*2.0;
		double t71 = t50*y0;
		double t73 = t23*t45*t49*(1.0/4.0);
		double t72 = t71-t73;
		double t74 = t45*t55*y0*(1.0/2.0);
		double t75 = t45*t49*y3*(1.0/2.0);
		double t76 = y0*y3*2.0;
		double t77 = y1*y3*2.0;
		double t78 = y4*y5*2.0;
		double t79 = t76+t77+t78;
		double t80 = t23*t45*t79*(1.0/4.0);
		double t107 = t23*t49*t55*t63*(3.0/8.0);
		double t81 = t74+t75+t80-t107;
		double t82 = t23*t28*t42*t81*(1.0/4.0);
		double t83 = t39*t42*(1.0/4.0);
		double t84 = t83-1.0;
		double t85 = y4*y4;
		double t86 = t85*2.0;
		double t87 = y5*y5;
		double t88 = t87*2.0;
		double t89 = t45*t61*y0*(1.0/2.0);
		double t90 = t45*t49*y4*(1.0/2.0);
		double t91 = y0*y4*2.0;
		double t92 = y2*y4*2.0;
		double t93 = y3*y5*2.0;
		double t94 = t91+t92+t93;
		double t95 = t23*t45*t94*(1.0/4.0);
		double t115 = t23*t49*t61*t63*(3.0/8.0);
		double t96 = t89+t90+t95-t115;
		double t97 = t23*t29*t42*t96*(1.0/4.0);
		double t98 = t45*t55*y4*(1.0/2.0);
		double t99 = t45*t61*y3*(1.0/2.0);
		double t100 = y1*y5*2.0;
		double t101 = y2*y5*2.0;
		double t102 = y3*y4*2.0;
		double t103 = t100+t101+t102;
		double t104 = t23*t45*t103*(1.0/4.0);
		double t124 = t23*t55*t61*t63*(3.0/8.0);
		double t105 = t98+t99+t104-t124;
		double t106 = t28*t29*t42*t105*(1.0/4.0);
		double t108 = t45*t49*y0;
		double t109 = y0*y0;
		double t110 = t109*2.0;
		double t111 = t70+t86+t110;
		double t112 = t23*t45*t111*(1.0/4.0);
		double t113 = t49*t49;
		double t116 = t23*t63*t113*(3.0/8.0);
		double t114 = t108+t112-t116;
		double t117 = t45*t55*y3;
		double t118 = y1*y1;
		double t119 = t118*2.0;
		double t120 = t70+t88+t119;
		double t121 = t23*t45*t120*(1.0/4.0);
		double t122 = t55*t55;
		double t125 = t23*t63*t122*(3.0/8.0);
		double t123 = t117+t121-t125;
		double t126 = t45*t61*y4;
		double t127 = y2*y2;
		double t128 = t127*2.0;
		double t129 = t86+t88+t128;
		double t130 = t23*t45*t129*(1.0/4.0);
		double t131 = t61*t61;
		double t133 = t23*t63*t131*(3.0/8.0);
		double t132 = t126+t130-t133;
		return -t44*(t82+t97+t44*t114+t72*(t27*t49*t65*(1.0/4.0)-t23*t42*y0)-t23*t42*t56*y3*(1.0/2.0)-t28*t42*t56*y0*(1.0/2.0)-t23*t42*t62*y4*(1.0/2.0)-t29*t42*t62*y0*(1.0/2.0)+t23*t28*t49*t56*t65*(1.0/4.0)+t23*t29*t49*t62*t65*(1.0/4.0))-t68*(t82+t106+t68*t123+t56*(t34*t55*t65*(1.0/4.0)-t28*t42*y1)-t29*t42*t62*y1*(1.0/2.0)-t28*t42*t62*y5*(1.0/2.0)-t23*t42*t72*y1*(1.0/2.0)-t28*t42*t72*y3*(1.0/2.0)+t28*t29*t55*t62*t65*(1.0/4.0)+t23*t28*t55*t65*t72*(1.0/4.0))-t84*(t97+t106+t84*t132+t62*(t39*t61*t65*(1.0/4.0)-t29*t42*y2)-t28*t42*t56*y2*(1.0/2.0)-t29*t42*t56*y5*(1.0/2.0)-t23*t42*t72*y2*(1.0/2.0)-t29*t42*t72*y4*(1.0/2.0)+t28*t29*t56*t61*t65*(1.0/4.0)+t23*t29*t61*t65*t72*(1.0/4.0))-t23*t28*t42*(t44*t81+t72*(t27*t55*t65*(1.0/4.0)-t23*t42*y3)+t23*t29*t42*t105*(1.0/4.0)+t23*t28*t42*t123*(1.0/4.0)-t23*t42*t56*y1*(1.0/2.0)-t28*t42*t56*y3*(1.0/2.0)-t23*t42*t62*y5*(1.0/2.0)-t29*t42*t62*y3*(1.0/2.0)+t23*t28*t55*t56*t65*(1.0/4.0)+t23*t29*t55*t62*t65*(1.0/4.0))*(1.0/4.0)-t23*t28*t42*(t68*t81+t56*(t34*t49*t65*(1.0/4.0)-t28*t42*y3)+t28*t29*t42*t96*(1.0/4.0)+t23*t28*t42*t114*(1.0/4.0)-t28*t42*t62*y4*(1.0/2.0)-t29*t42*t62*y3*(1.0/2.0)-t23*t42*t72*y3*(1.0/2.0)-t28*t42*t72*y0*(1.0/2.0)+t28*t29*t49*t62*t65*(1.0/4.0)+t23*t28*t49*t65*t72*(1.0/4.0))*(1.0/4.0)-t23*t29*t42*(t44*t96+t72*(t27*t61*t65*(1.0/4.0)-t23*t42*y4)+t23*t28*t42*t105*(1.0/4.0)+t23*t29*t42*t132*(1.0/4.0)-t23*t42*t56*y5*(1.0/2.0)-t23*t42*t62*y2*(1.0/2.0)-t28*t42*t56*y4*(1.0/2.0)-t29*t42*t62*y4*(1.0/2.0)+t23*t28*t56*t61*t65*(1.0/4.0)+t23*t29*t61*t62*t65*(1.0/4.0))*(1.0/4.0)-t23*t29*t42*(t84*t96+t62*(t39*t49*t65*(1.0/4.0)-t29*t42*y4)+t28*t29*t42*t81*(1.0/4.0)+t23*t29*t42*t114*(1.0/4.0)-t28*t42*t56*y4*(1.0/2.0)-t29*t42*t56*y3*(1.0/2.0)-t23*t42*t72*y4*(1.0/2.0)-t29*t42*t72*y0*(1.0/2.0)+t28*t29*t49*t56*t65*(1.0/4.0)+t23*t29*t49*t65*t72*(1.0/4.0))*(1.0/4.0)-t28*t29*t42*(t84*t105+t62*(t39*t55*t65*(1.0/4.0)-t29*t42*y5)+t23*t29*t42*t81*(1.0/4.0)+t28*t29*t42*t123*(1.0/4.0)-t29*t42*t56*y1*(1.0/2.0)-t28*t42*t56*y5*(1.0/2.0)-t23*t42*t72*y5*(1.0/2.0)-t29*t42*t72*y3*(1.0/2.0)+t28*t29*t55*t56*t65*(1.0/4.0)+t23*t29*t55*t65*t72*(1.0/4.0))*(1.0/4.0)-t28*t29*t42*(t68*t105+t56*(t34*t61*t65*(1.0/4.0)-t28*t42*y5)+t23*t28*t42*t96*(1.0/4.0)+t28*t29*t42*t132*(1.0/4.0)-t28*t42*t62*y2*(1.0/2.0)-t29*t42*t62*y5*(1.0/2.0)-t23*t42*t72*y5*(1.0/2.0)-t28*t42*t72*y4*(1.0/2.0)+t28*t29*t61*t62*t65*(1.0/4.0)+t23*t28*t61*t65*t72*(1.0/4.0))*(1.0/4.0);
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
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t136 = x0*y0*2.0;
		double t137 = x1*y3*2.0;
		double t138 = x2*y4*2.0;
		double t135 = t136+t137+t138;
		double t139 = t135*t135;
		double t143 = x1*y1*2.0;
		double t144 = x0*y3*2.0;
		double t145 = x2*y5*2.0;
		double t140 = t143+t144+t145;
		double t148 = x0*y4*2.0;
		double t149 = x2*y2*2.0;
		double t150 = x1*y5*2.0;
		double t141 = t148+t149+t150;
		double t142 = t139*(1.0/4.0);
		double t146 = t140*t140;
		double t147 = t146*(1.0/4.0);
		double t151 = t141*t141;
		double t152 = t151*(1.0/4.0);
		double t153 = t142+t147+t152;
		double t154 = 1.0/t153;
		double t155 = t139*t154*(1.0/4.0);
		double t156 = t155-1.0;
		double t157 = 1.0/pow(t153,3.0/2.0);
		double t158 = t135*y0;
		double t159 = t140*y3;
		double t160 = t141*y4;
		double t161 = t158+t159+t160;
		double t162 = 1.0/sqrt(t153);
		double t163 = t162*y1;
		double t164 = t140*y1;
		double t165 = t135*y3;
		double t166 = t141*y5;
		double t167 = t164+t165+t166;
		double t176 = t140*t157*t167*(1.0/4.0);
		double t168 = t163-t176;
		double t169 = t162*y5;
		double t170 = t135*y4;
		double t171 = t141*y2;
		double t172 = t140*y5;
		double t173 = t170+t171+t172;
		double t178 = t140*t157*t173*(1.0/4.0);
		double t174 = t169-t178;
		double t175 = 1.0/pow(t153,5.0/2.0);
		double t177 = 1.0/(t153*t153);
		double t179 = t146*t154*(1.0/4.0);
		double t180 = t179-1.0;
		double t181 = y3*y3;
		double t182 = t181*2.0;
		double t183 = t162*y3;
		double t185 = t140*t157*t161*(1.0/4.0);
		double t184 = t183-t185;
		double t186 = t157*t161*y1*(1.0/2.0);
		double t187 = t157*t167*y3*(1.0/2.0);
		double t188 = y0*y3*2.0;
		double t189 = y1*y3*2.0;
		double t190 = y4*y5*2.0;
		double t191 = t188+t189+t190;
		double t192 = t140*t157*t191*(1.0/4.0);
		double t219 = t140*t161*t167*t175*(3.0/8.0);
		double t193 = t186+t187+t192-t219;
		double t194 = t135*t140*t154*t193*(1.0/4.0);
		double t195 = t151*t154*(1.0/4.0);
		double t196 = t195-1.0;
		double t197 = y4*y4;
		double t198 = t197*2.0;
		double t199 = y5*y5;
		double t200 = t199*2.0;
		double t201 = t157*t161*y5*(1.0/2.0);
		double t202 = t157*t173*y3*(1.0/2.0);
		double t203 = y0*y4*2.0;
		double t204 = y2*y4*2.0;
		double t205 = y3*y5*2.0;
		double t206 = t203+t204+t205;
		double t207 = t140*t157*t206*(1.0/4.0);
		double t227 = t140*t161*t173*t175*(3.0/8.0);
		double t208 = t201+t202+t207-t227;
		double t209 = t135*t141*t154*t208*(1.0/4.0);
		double t210 = t157*t173*y1*(1.0/2.0);
		double t211 = t157*t167*y5*(1.0/2.0);
		double t212 = y1*y5*2.0;
		double t213 = y2*y5*2.0;
		double t214 = y3*y4*2.0;
		double t215 = t212+t213+t214;
		double t216 = t140*t157*t215*(1.0/4.0);
		double t235 = t140*t167*t173*t175*(3.0/8.0);
		double t217 = t210+t211+t216-t235;
		double t218 = t140*t141*t154*t217*(1.0/4.0);
		double t220 = t157*t161*y3;
		double t221 = y0*y0;
		double t222 = t221*2.0;
		double t223 = t182+t198+t222;
		double t224 = t140*t157*t223*(1.0/4.0);
		double t225 = t161*t161;
		double t236 = t140*t175*t225*(3.0/8.0);
		double t226 = t220+t224-t236;
		double t228 = t157*t167*y1;
		double t229 = y1*y1;
		double t230 = t229*2.0;
		double t231 = t182+t200+t230;
		double t232 = t140*t157*t231*(1.0/4.0);
		double t233 = t167*t167;
		double t237 = t140*t175*t233*(3.0/8.0);
		double t234 = t228+t232-t237;
		double t238 = t157*t173*y5;
		double t239 = y2*y2;
		double t240 = t239*2.0;
		double t241 = t198+t200+t240;
		double t242 = t140*t157*t241*(1.0/4.0);
		double t243 = t173*t173;
		double t245 = t140*t175*t243*(3.0/8.0);
		double t244 = t238+t242-t245;
		return -t156*(t194+t209+t156*t226+t184*(t139*t161*t177*(1.0/4.0)-t135*t154*y0)-t135*t154*t168*y3*(1.0/2.0)-t140*t154*t168*y0*(1.0/2.0)-t135*t154*t174*y4*(1.0/2.0)-t141*t154*t174*y0*(1.0/2.0)+t135*t140*t161*t168*t177*(1.0/4.0)+t135*t141*t161*t174*t177*(1.0/4.0))-t180*(t194+t218+t180*t234+t168*(t146*t167*t177*(1.0/4.0)-t140*t154*y1)-t141*t154*t174*y1*(1.0/2.0)-t140*t154*t174*y5*(1.0/2.0)-t135*t154*t184*y1*(1.0/2.0)-t140*t154*t184*y3*(1.0/2.0)+t140*t141*t167*t174*t177*(1.0/4.0)+t135*t140*t167*t177*t184*(1.0/4.0))-t196*(t209+t218+t196*t244+t174*(t151*t173*t177*(1.0/4.0)-t141*t154*y2)-t140*t154*t168*y2*(1.0/2.0)-t141*t154*t168*y5*(1.0/2.0)-t135*t154*t184*y2*(1.0/2.0)-t141*t154*t184*y4*(1.0/2.0)+t140*t141*t168*t173*t177*(1.0/4.0)+t135*t141*t173*t177*t184*(1.0/4.0))-t135*t140*t154*(t156*t193+t184*(t139*t167*t177*(1.0/4.0)-t135*t154*y3)+t135*t141*t154*t217*(1.0/4.0)+t135*t140*t154*t234*(1.0/4.0)-t135*t154*t168*y1*(1.0/2.0)-t140*t154*t168*y3*(1.0/2.0)-t135*t154*t174*y5*(1.0/2.0)-t141*t154*t174*y3*(1.0/2.0)+t135*t140*t167*t168*t177*(1.0/4.0)+t135*t141*t167*t174*t177*(1.0/4.0))*(1.0/4.0)-t135*t140*t154*(t180*t193+t168*(t146*t161*t177*(1.0/4.0)-t140*t154*y3)+t140*t141*t154*t208*(1.0/4.0)+t135*t140*t154*t226*(1.0/4.0)-t140*t154*t174*y4*(1.0/2.0)-t141*t154*t174*y3*(1.0/2.0)-t135*t154*t184*y3*(1.0/2.0)-t140*t154*t184*y0*(1.0/2.0)+t140*t141*t161*t174*t177*(1.0/4.0)+t135*t140*t161*t177*t184*(1.0/4.0))*(1.0/4.0)-t135*t141*t154*(t156*t208+t184*(t139*t173*t177*(1.0/4.0)-t135*t154*y4)+t135*t140*t154*t217*(1.0/4.0)+t135*t141*t154*t244*(1.0/4.0)-t135*t154*t168*y5*(1.0/2.0)-t135*t154*t174*y2*(1.0/2.0)-t140*t154*t168*y4*(1.0/2.0)-t141*t154*t174*y4*(1.0/2.0)+t135*t140*t168*t173*t177*(1.0/4.0)+t135*t141*t173*t174*t177*(1.0/4.0))*(1.0/4.0)-t135*t141*t154*(t196*t208+t174*(t151*t161*t177*(1.0/4.0)-t141*t154*y4)+t140*t141*t154*t193*(1.0/4.0)+t135*t141*t154*t226*(1.0/4.0)-t140*t154*t168*y4*(1.0/2.0)-t141*t154*t168*y3*(1.0/2.0)-t135*t154*t184*y4*(1.0/2.0)-t141*t154*t184*y0*(1.0/2.0)+t140*t141*t161*t168*t177*(1.0/4.0)+t135*t141*t161*t177*t184*(1.0/4.0))*(1.0/4.0)-t140*t141*t154*(t196*t217+t174*(t151*t167*t177*(1.0/4.0)-t141*t154*y5)+t135*t141*t154*t193*(1.0/4.0)+t140*t141*t154*t234*(1.0/4.0)-t141*t154*t168*y1*(1.0/2.0)-t140*t154*t168*y5*(1.0/2.0)-t135*t154*t184*y5*(1.0/2.0)-t141*t154*t184*y3*(1.0/2.0)+t140*t141*t167*t168*t177*(1.0/4.0)+t135*t141*t167*t177*t184*(1.0/4.0))*(1.0/4.0)-t140*t141*t154*(t180*t217+t168*(t146*t173*t177*(1.0/4.0)-t140*t154*y5)+t135*t140*t154*t208*(1.0/4.0)+t140*t141*t154*t244*(1.0/4.0)-t140*t154*t174*y2*(1.0/2.0)-t141*t154*t174*y5*(1.0/2.0)-t135*t154*t184*y5*(1.0/2.0)-t140*t154*t184*y4*(1.0/2.0)+t140*t141*t173*t174*t177*(1.0/4.0)+t135*t140*t173*t177*t184*(1.0/4.0))*(1.0/4.0);
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
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t248 = x0*y0*2.0;
		double t249 = x1*y3*2.0;
		double t250 = x2*y4*2.0;
		double t247 = t248+t249+t250;
		double t251 = t247*t247;
		double t255 = x1*y1*2.0;
		double t256 = x0*y3*2.0;
		double t257 = x2*y5*2.0;
		double t252 = t255+t256+t257;
		double t260 = x0*y4*2.0;
		double t261 = x2*y2*2.0;
		double t262 = x1*y5*2.0;
		double t253 = t260+t261+t262;
		double t254 = t251*(1.0/4.0);
		double t258 = t252*t252;
		double t259 = t258*(1.0/4.0);
		double t263 = t253*t253;
		double t264 = t263*(1.0/4.0);
		double t265 = t254+t259+t264;
		double t266 = 1.0/t265;
		double t267 = t251*t266*(1.0/4.0);
		double t268 = t267-1.0;
		double t269 = 1.0/pow(t265,3.0/2.0);
		double t270 = t247*y0;
		double t271 = t252*y3;
		double t272 = t253*y4;
		double t273 = t270+t271+t272;
		double t274 = 1.0/sqrt(t265);
		double t275 = t274*y5;
		double t276 = t252*y1;
		double t277 = t247*y3;
		double t278 = t253*y5;
		double t279 = t276+t277+t278;
		double t288 = t253*t269*t279*(1.0/4.0);
		double t280 = t275-t288;
		double t281 = t274*y2;
		double t282 = t247*y4;
		double t283 = t253*y2;
		double t284 = t252*y5;
		double t285 = t282+t283+t284;
		double t290 = t253*t269*t285*(1.0/4.0);
		double t286 = t281-t290;
		double t287 = 1.0/pow(t265,5.0/2.0);
		double t289 = 1.0/(t265*t265);
		double t291 = t258*t266*(1.0/4.0);
		double t292 = t291-1.0;
		double t293 = y3*y3;
		double t294 = t293*2.0;
		double t295 = t274*y4;
		double t297 = t253*t269*t273*(1.0/4.0);
		double t296 = t295-t297;
		double t298 = t269*t273*y5*(1.0/2.0);
		double t299 = t269*t279*y4*(1.0/2.0);
		double t300 = y0*y3*2.0;
		double t301 = y1*y3*2.0;
		double t302 = y4*y5*2.0;
		double t303 = t300+t301+t302;
		double t304 = t253*t269*t303*(1.0/4.0);
		double t331 = t253*t273*t279*t287*(3.0/8.0);
		double t305 = t298+t299+t304-t331;
		double t306 = t247*t252*t266*t305*(1.0/4.0);
		double t307 = t263*t266*(1.0/4.0);
		double t308 = t307-1.0;
		double t309 = y4*y4;
		double t310 = t309*2.0;
		double t311 = y5*y5;
		double t312 = t311*2.0;
		double t313 = t269*t273*y2*(1.0/2.0);
		double t314 = t269*t285*y4*(1.0/2.0);
		double t315 = y0*y4*2.0;
		double t316 = y2*y4*2.0;
		double t317 = y3*y5*2.0;
		double t318 = t315+t316+t317;
		double t319 = t253*t269*t318*(1.0/4.0);
		double t339 = t253*t273*t285*t287*(3.0/8.0);
		double t320 = t313+t314+t319-t339;
		double t321 = t247*t253*t266*t320*(1.0/4.0);
		double t322 = t269*t279*y2*(1.0/2.0);
		double t323 = t269*t285*y5*(1.0/2.0);
		double t324 = y1*y5*2.0;
		double t325 = y2*y5*2.0;
		double t326 = y3*y4*2.0;
		double t327 = t324+t325+t326;
		double t328 = t253*t269*t327*(1.0/4.0);
		double t347 = t253*t279*t285*t287*(3.0/8.0);
		double t329 = t322+t323+t328-t347;
		double t330 = t252*t253*t266*t329*(1.0/4.0);
		double t332 = t269*t273*y4;
		double t333 = y0*y0;
		double t334 = t333*2.0;
		double t335 = t294+t310+t334;
		double t336 = t253*t269*t335*(1.0/4.0);
		double t337 = t273*t273;
		double t348 = t253*t287*t337*(3.0/8.0);
		double t338 = t332+t336-t348;
		double t340 = t269*t279*y5;
		double t341 = y1*y1;
		double t342 = t341*2.0;
		double t343 = t294+t312+t342;
		double t344 = t253*t269*t343*(1.0/4.0);
		double t345 = t279*t279;
		double t357 = t253*t287*t345*(3.0/8.0);
		double t346 = t340+t344-t357;
		double t349 = t269*t285*y2;
		double t350 = y2*y2;
		double t351 = t350*2.0;
		double t352 = t310+t312+t351;
		double t353 = t253*t269*t352*(1.0/4.0);
		double t354 = t285*t285;
		double t356 = t253*t287*t354*(3.0/8.0);
		double t355 = t349+t353-t356;
		return -t268*(t306+t321+t268*t338+t296*(t251*t273*t289*(1.0/4.0)-t247*t266*y0)-t247*t266*t280*y3*(1.0/2.0)-t252*t266*t280*y0*(1.0/2.0)-t247*t266*t286*y4*(1.0/2.0)-t253*t266*t286*y0*(1.0/2.0)+t247*t252*t273*t280*t289*(1.0/4.0)+t247*t253*t273*t286*t289*(1.0/4.0))-t292*(t306+t330+t292*t346+t280*(t258*t279*t289*(1.0/4.0)-t252*t266*y1)-t253*t266*t286*y1*(1.0/2.0)-t252*t266*t286*y5*(1.0/2.0)-t247*t266*t296*y1*(1.0/2.0)-t252*t266*t296*y3*(1.0/2.0)+t252*t253*t279*t286*t289*(1.0/4.0)+t247*t252*t279*t289*t296*(1.0/4.0))-t308*(t321+t330+t308*t355+t286*(t263*t285*t289*(1.0/4.0)-t253*t266*y2)-t252*t266*t280*y2*(1.0/2.0)-t253*t266*t280*y5*(1.0/2.0)-t247*t266*t296*y2*(1.0/2.0)-t253*t266*t296*y4*(1.0/2.0)+t252*t253*t280*t285*t289*(1.0/4.0)+t247*t253*t285*t289*t296*(1.0/4.0))-t247*t252*t266*(t268*t305+t296*(t251*t279*t289*(1.0/4.0)-t247*t266*y3)+t247*t253*t266*t329*(1.0/4.0)+t247*t252*t266*t346*(1.0/4.0)-t247*t266*t280*y1*(1.0/2.0)-t252*t266*t280*y3*(1.0/2.0)-t247*t266*t286*y5*(1.0/2.0)-t253*t266*t286*y3*(1.0/2.0)+t247*t252*t279*t280*t289*(1.0/4.0)+t247*t253*t279*t286*t289*(1.0/4.0))*(1.0/4.0)-t247*t252*t266*(t292*t305+t280*(t258*t273*t289*(1.0/4.0)-t252*t266*y3)+t252*t253*t266*t320*(1.0/4.0)+t247*t252*t266*t338*(1.0/4.0)-t252*t266*t286*y4*(1.0/2.0)-t253*t266*t286*y3*(1.0/2.0)-t247*t266*t296*y3*(1.0/2.0)-t252*t266*t296*y0*(1.0/2.0)+t252*t253*t273*t286*t289*(1.0/4.0)+t247*t252*t273*t289*t296*(1.0/4.0))*(1.0/4.0)-t247*t253*t266*(t268*t320+t296*(t251*t285*t289*(1.0/4.0)-t247*t266*y4)+t247*t252*t266*t329*(1.0/4.0)+t247*t253*t266*t355*(1.0/4.0)-t247*t266*t280*y5*(1.0/2.0)-t247*t266*t286*y2*(1.0/2.0)-t252*t266*t280*y4*(1.0/2.0)-t253*t266*t286*y4*(1.0/2.0)+t247*t252*t280*t285*t289*(1.0/4.0)+t247*t253*t285*t286*t289*(1.0/4.0))*(1.0/4.0)-t247*t253*t266*(t308*t320+t286*(t263*t273*t289*(1.0/4.0)-t253*t266*y4)+t252*t253*t266*t305*(1.0/4.0)+t247*t253*t266*t338*(1.0/4.0)-t252*t266*t280*y4*(1.0/2.0)-t253*t266*t280*y3*(1.0/2.0)-t247*t266*t296*y4*(1.0/2.0)-t253*t266*t296*y0*(1.0/2.0)+t252*t253*t273*t280*t289*(1.0/4.0)+t247*t253*t273*t289*t296*(1.0/4.0))*(1.0/4.0)-t252*t253*t266*(t308*t329+t286*(t263*t279*t289*(1.0/4.0)-t253*t266*y5)+t247*t253*t266*t305*(1.0/4.0)+t252*t253*t266*t346*(1.0/4.0)-t253*t266*t280*y1*(1.0/2.0)-t252*t266*t280*y5*(1.0/2.0)-t247*t266*t296*y5*(1.0/2.0)-t253*t266*t296*y3*(1.0/2.0)+t252*t253*t279*t280*t289*(1.0/4.0)+t247*t253*t279*t289*t296*(1.0/4.0))*(1.0/4.0)-t252*t253*t266*(t292*t329+t280*(t258*t285*t289*(1.0/4.0)-t252*t266*y5)+t247*t252*t266*t320*(1.0/4.0)+t252*t253*t266*t355*(1.0/4.0)-t252*t266*t286*y2*(1.0/2.0)-t253*t266*t286*y5*(1.0/2.0)-t247*t266*t296*y5*(1.0/2.0)-t252*t266*t296*y4*(1.0/2.0)+t252*t253*t285*t286*t289*(1.0/4.0)+t247*t252*t285*t289*t296*(1.0/4.0))*(1.0/4.0);
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
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t605 = x0*y0*2.0;
		double t606 = x1*y3*2.0;
		double t607 = x2*y4*2.0;
		double t604 = t605+t606+t607;
		double t608 = t604*t604;
		double t612 = x1*y1*2.0;
		double t613 = x0*y3*2.0;
		double t614 = x2*y5*2.0;
		double t609 = t612+t613+t614;
		double t617 = x0*y4*2.0;
		double t618 = x2*y2*2.0;
		double t619 = x1*y5*2.0;
		double t610 = t617+t618+t619;
		double t611 = t608*(1.0/4.0);
		double t615 = t609*t609;
		double t616 = t615*(1.0/4.0);
		double t620 = t610*t610;
		double t621 = t620*(1.0/4.0);
		double t622 = t611+t616+t621;
		double t623 = 1.0/t622;
		double t624 = 1.0/sqrt(t622);
		double t625 = 1.0/pow(t622,3.0/2.0);
		double t626 = t609*y1;
		double t627 = t604*y3;
		double t628 = t610*y5;
		double t629 = t626+t627+t628;
		double t630 = t624*y3;
		double t631 = t604*y0;
		double t632 = t609*y3;
		double t633 = t610*y4;
		double t634 = t631+t632+t633;
		double t635 = t604*y4;
		double t636 = t610*y2;
		double t637 = t609*y5;
		double t638 = t635+t636+t637;
		double t639 = t624*y4;
		double t640 = t624*y5;
		double t642 = t608*t623*(1.0/4.0);
		double t643 = t642-1.0;
		double t644 = t624*y0;
		double t645 = t604*t625*t634*(1.0/4.0);
		double t646 = t644-t645;
		double t647 = t643*t646;
		double t648 = t604*t625*t629*(1.0/4.0);
		double t649 = t630-t648;
		double t650 = t604*t609*t623*t649*(1.0/4.0);
		double t651 = t604*t625*t638*(1.0/4.0);
		double t652 = t639-t651;
		double t653 = t604*t610*t623*t652*(1.0/4.0);
		double t655 = t615*t623*(1.0/4.0);
		double t656 = t655-1.0;
		double t658 = t609*t625*t634*(1.0/4.0);
		double t659 = t630-t658;
		double t660 = t624*y1;
		double t661 = t609*t625*t629*(1.0/4.0);
		double t662 = t660-t661;
		double t663 = t609*t625*t638*(1.0/4.0);
		double t664 = t640-t663;
		double t666 = t620*t623*(1.0/4.0);
		double t667 = t666-1.0;
		double t669 = t656*t662;
		double t670 = t604*t609*t623*t659*(1.0/4.0);
		double t671 = t609*t610*t623*t664*(1.0/4.0);
		double t673 = t610*t625*t634*(1.0/4.0);
		double t674 = t639-t673;
		double t675 = t610*t625*t629*(1.0/4.0);
		double t676 = t640-t675;
		double t677 = t624*y2;
		double t678 = t610*t625*t638*(1.0/4.0);
		double t679 = t677-t678;
		double t683 = t667*t679;
		double t684 = t604*t610*t623*t674*(1.0/4.0);
		double t685 = t609*t610*t623*t676*(1.0/4.0);
		double t641 = t647+t650+t653+t669+t670+t671+t683+t684+t685;
		double t654 = t647+t650+t653;
		double t657 = t649*t656+t604*t609*t623*t646*(1.0/4.0)+t609*t610*t623*t652*(1.0/4.0);
		double t665 = t643*t659+t604*t609*t623*t662*(1.0/4.0)+t604*t610*t623*t664*(1.0/4.0);
		double t668 = t652*t667+t604*t610*t623*t646*(1.0/4.0)+t609*t610*t623*t649*(1.0/4.0);
		double t672 = t669+t670+t671;
		double t680 = t643*t674+t604*t609*t623*t676*(1.0/4.0)+t604*t610*t623*t679*(1.0/4.0);
		double t681 = t664*t667+t604*t610*t623*t659*(1.0/4.0)+t609*t610*t623*t662*(1.0/4.0);
		double t682 = t656*t676+t604*t609*t623*t674*(1.0/4.0)+t609*t610*t623*t679*(1.0/4.0);
		double t686 = t683+t684+t685;
		return (t641*t641)*(1.0/2.0)-(t654*t654)*(1.0/2.0)-(t657*t657)*(1.0/2.0)-(t665*t665)*(1.0/2.0)-(t668*t668)*(1.0/2.0)-(t672*t672)*(1.0/2.0)-(t680*t680)*(1.0/2.0)-(t681*t681)*(1.0/2.0)-(t682*t682)*(1.0/2.0)-(t686*t686)*(1.0/2.0);
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
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		double t567 = x0*y0*2.0;
		double t568 = x1*y3*2.0;
		double t569 = x2*y4*2.0;
		double t566 = t567+t568+t569;
		double t570 = t566*t566;
		double t574 = x1*y1*2.0;
		double t575 = x0*y3*2.0;
		double t576 = x2*y5*2.0;
		double t571 = t574+t575+t576;
		double t579 = x0*y4*2.0;
		double t580 = x2*y2*2.0;
		double t581 = x1*y5*2.0;
		double t572 = t579+t580+t581;
		double t573 = t570*(1.0/4.0);
		double t577 = t571*t571;
		double t578 = t577*(1.0/4.0);
		double t582 = t572*t572;
		double t583 = t582*(1.0/4.0);
		double t584 = t573+t578+t583;
		double t585 = 1.0/t584;
		double t586 = 1.0/sqrt(t584);
		double t587 = 1.0/pow(t584,3.0/2.0);
		double t588 = t571*y1;
		double t589 = t566*y3;
		double t590 = t572*y5;
		double t591 = t588+t589+t590;
		double t592 = t586*y3;
		double t593 = t566*y0;
		double t594 = t571*y3;
		double t595 = t572*y4;
		double t596 = t593+t594+t595;
		double t597 = t566*y4;
		double t598 = t572*y2;
		double t599 = t571*y5;
		double t600 = t597+t598+t599;
		double t601 = t586*y4;
		double t602 = t586*y5;
		return -(t586*y0-t566*t587*t596*(1.0/4.0))*(t570*t585*(1.0/4.0)-1.0)-(t586*y1-t571*t587*t591*(1.0/4.0))*(t577*t585*(1.0/4.0)-1.0)-(t586*y2-t572*t587*t600*(1.0/4.0))*(t582*t585*(1.0/4.0)-1.0)-t566*t571*t585*(t592-t566*t587*t591*(1.0/4.0))*(1.0/4.0)-t566*t571*t585*(t592-t571*t587*t596*(1.0/4.0))*(1.0/4.0)-t566*t572*t585*(t601-t566*t587*t600*(1.0/4.0))*(1.0/4.0)-t566*t572*t585*(t601-t572*t587*t596*(1.0/4.0))*(1.0/4.0)-t571*t572*t585*(t602-t572*t587*t591*(1.0/4.0))*(1.0/4.0)-t571*t572*t585*(t602-t571*t587*t600*(1.0/4.0))*(1.0/4.0);
	}

	static inline WorldVector<double> getGradN0(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		WorldVector<double> res;
		// gradN0[0]
		{
			double t360 = x0*y0*2.0;
			double t361 = x1*y3*2.0;
			double t362 = x2*y4*2.0;
			double t359 = t360+t361+t362;
			double t363 = t359*t359;
			double t367 = x1*y1*2.0;
			double t368 = x0*y3*2.0;
			double t369 = x2*y5*2.0;
			double t364 = t367+t368+t369;
			double t372 = x0*y4*2.0;
			double t373 = x2*y2*2.0;
			double t374 = x1*y5*2.0;
			double t365 = t372+t373+t374;
			double t366 = t363*(1.0/4.0);
			double t370 = t364*t364;
			double t371 = t370*(1.0/4.0);
			double t375 = t365*t365;
			double t376 = t375*(1.0/4.0);
			double t377 = t366+t371+t376;
			double t378 = 1.0/sqrt(t377);
			double t379 = 1.0/pow(t377,3.0/2.0);
			double t380 = 1.0/t377;
			res[0] = -(t378*y0-t359*t379*(t359*y0+t364*y3+t365*y4)*(1.0/4.0))*(t363*t380*(1.0/4.0)-1.0)-t359*t364*t380*(t378*y3-t359*t379*(t359*y3+t364*y1+t365*y5)*(1.0/4.0))*(1.0/4.0)-t359*t365*t380*(t378*y4-t359*t379*(t359*y4+t365*y2+t364*y5)*(1.0/4.0))*(1.0/4.0);
		}
		// gradN0[1]
		{
			double t384 = x1*y1*2.0;
			double t385 = x0*y3*2.0;
			double t386 = x2*y5*2.0;
			double t382 = t384+t385+t386;
			double t389 = x0*y0*2.0;
			double t390 = x1*y3*2.0;
			double t391 = x2*y4*2.0;
			double t383 = t389+t390+t391;
			double t387 = t382*t382;
			double t395 = x0*y4*2.0;
			double t396 = x2*y2*2.0;
			double t397 = x1*y5*2.0;
			double t388 = t395+t396+t397;
			double t392 = t383*t383;
			double t393 = t392*(1.0/4.0);
			double t394 = t387*(1.0/4.0);
			double t398 = t388*t388;
			double t399 = t398*(1.0/4.0);
			double t400 = t393+t394+t399;
			double t401 = 1.0/sqrt(t400);
			double t402 = 1.0/pow(t400,3.0/2.0);
			double t403 = 1.0/t400;
			res[1] = -(t401*y3-t383*t402*(t382*y1+t383*y3+t388*y5)*(1.0/4.0))*(t387*t403*(1.0/4.0)-1.0)-t382*t383*t403*(t401*y0-t383*t402*(t383*y0+t382*y3+t388*y4)*(1.0/4.0))*(1.0/4.0)-t382*t388*t403*(t401*y4-t383*t402*(t382*y5+t383*y4+t388*y2)*(1.0/4.0))*(1.0/4.0);
		}
		// gradN0[2]
		{
			double t408 = x0*y4*2.0;
			double t409 = x2*y2*2.0;
			double t410 = x1*y5*2.0;
			double t405 = t408+t409+t410;
			double t412 = x0*y0*2.0;
			double t413 = x1*y3*2.0;
			double t414 = x2*y4*2.0;
			double t406 = t412+t413+t414;
			double t417 = x1*y1*2.0;
			double t418 = x0*y3*2.0;
			double t419 = x2*y5*2.0;
			double t407 = t417+t418+t419;
			double t411 = t405*t405;
			double t415 = t406*t406;
			double t416 = t415*(1.0/4.0);
			double t420 = t407*t407;
			double t421 = t420*(1.0/4.0);
			double t422 = t411*(1.0/4.0);
			double t423 = t416+t421+t422;
			double t424 = 1.0/sqrt(t423);
			double t425 = 1.0/pow(t423,3.0/2.0);
			double t426 = 1.0/t423;
			res[2] = -(t424*y4-t406*t425*(t405*y2+t406*y4+t407*y5)*(1.0/4.0))*(t411*t426*(1.0/4.0)-1.0)-t405*t406*t426*(t424*y0-t406*t425*(t406*y0+t405*y4+t407*y3)*(1.0/4.0))*(1.0/4.0)-t405*t407*t426*(t424*y3-t406*t425*(t407*y1+t406*y3+t405*y5)*(1.0/4.0))*(1.0/4.0);
		}
		return res ;
	}

	static inline WorldVector<double> getGradN1(
			const mtl::dense_vector<double> y ,
			const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		WorldVector<double> res;
		// gradN1[0]
		{
			double t429 = x0*y0*2.0;
			double t430 = x1*y3*2.0;
			double t431 = x2*y4*2.0;
			double t428 = t429+t430+t431;
			double t432 = t428*t428;
			double t436 = x1*y1*2.0;
			double t437 = x0*y3*2.0;
			double t438 = x2*y5*2.0;
			double t433 = t436+t437+t438;
			double t441 = x0*y4*2.0;
			double t442 = x2*y2*2.0;
			double t443 = x1*y5*2.0;
			double t434 = t441+t442+t443;
			double t435 = t432*(1.0/4.0);
			double t439 = t433*t433;
			double t440 = t439*(1.0/4.0);
			double t444 = t434*t434;
			double t445 = t444*(1.0/4.0);
			double t446 = t435+t440+t445;
			double t447 = 1.0/sqrt(t446);
			double t448 = 1.0/pow(t446,3.0/2.0);
			double t449 = 1.0/t446;
			res[0] = -(t447*y3-t433*t448*(t428*y0+t433*y3+t434*y4)*(1.0/4.0))*(t432*t449*(1.0/4.0)-1.0)-t428*t433*t449*(t447*y1-t433*t448*(t428*y3+t433*y1+t434*y5)*(1.0/4.0))*(1.0/4.0)-t428*t434*t449*(t447*y5-t433*t448*(t428*y4+t434*y2+t433*y5)*(1.0/4.0))*(1.0/4.0);
		}
		// gradN1[1]
		{
			double t453 = x1*y1*2.0;
			double t454 = x0*y3*2.0;
			double t455 = x2*y5*2.0;
			double t451 = t453+t454+t455;
			double t458 = x0*y0*2.0;
			double t459 = x1*y3*2.0;
			double t460 = x2*y4*2.0;
			double t452 = t458+t459+t460;
			double t456 = t451*t451;
			double t464 = x0*y4*2.0;
			double t465 = x2*y2*2.0;
			double t466 = x1*y5*2.0;
			double t457 = t464+t465+t466;
			double t461 = t452*t452;
			double t462 = t461*(1.0/4.0);
			double t463 = t456*(1.0/4.0);
			double t467 = t457*t457;
			double t468 = t467*(1.0/4.0);
			double t469 = t462+t463+t468;
			double t470 = 1.0/sqrt(t469);
			double t471 = 1.0/pow(t469,3.0/2.0);
			double t472 = 1.0/t469;
			res[1] = -(t470*y1-t451*t471*(t451*y1+t452*y3+t457*y5)*(1.0/4.0))*(t456*t472*(1.0/4.0)-1.0)-t451*t452*t472*(t470*y3-t451*t471*(t452*y0+t451*y3+t457*y4)*(1.0/4.0))*(1.0/4.0)-t451*t457*t472*(t470*y5-t451*t471*(t451*y5+t452*y4+t457*y2)*(1.0/4.0))*(1.0/4.0);
		}
		// gradN1[2]
		{
			double t477 = x0*y4*2.0;
			double t478 = x2*y2*2.0;
			double t479 = x1*y5*2.0;
			double t474 = t477+t478+t479;
			double t481 = x0*y0*2.0;
			double t482 = x1*y3*2.0;
			double t483 = x2*y4*2.0;
			double t475 = t481+t482+t483;
			double t486 = x1*y1*2.0;
			double t487 = x0*y3*2.0;
			double t488 = x2*y5*2.0;
			double t476 = t486+t487+t488;
			double t480 = t474*t474;
			double t484 = t475*t475;
			double t485 = t484*(1.0/4.0);
			double t489 = t476*t476;
			double t490 = t489*(1.0/4.0);
			double t491 = t480*(1.0/4.0);
			double t492 = t485+t490+t491;
			double t493 = 1.0/sqrt(t492);
			double t494 = 1.0/pow(t492,3.0/2.0);
			double t495 = 1.0/t492;
			res[2] = -(t493*y5-t476*t494*(t474*y2+t475*y4+t476*y5)*(1.0/4.0))*(t480*t495*(1.0/4.0)-1.0)-t474*t475*t495*(t493*y3-t476*t494*(t475*y0+t474*y4+t476*y3)*(1.0/4.0))*(1.0/4.0)-t474*t476*t495*(t493*y1-t476*t494*(t476*y1+t475*y3+t474*y5)*(1.0/4.0))*(1.0/4.0);
		}
		return res ;
	}

	static inline WorldVector<double> getGradN2(
				const mtl::dense_vector<double> y ,
				const WorldVector<double> &x )
	{
		double y0 = y[0] ;
		double y1 = y[1] ;
		double y2 = y[2] ;
		double y3 = y[3] ;
		double y4 = y[4] ;
		double y5 = y[5] ;
		double x0 = x[0] ;
		double x1 = x[1] ;
		double x2 = x[2] ;
		WorldVector<double> res;
		// gradN2[0]
		{
			double t498 = x0*y0*2.0;
			double t499 = x1*y3*2.0;
			double t500 = x2*y4*2.0;
			double t497 = t498+t499+t500;
			double t501 = t497*t497;
			double t505 = x1*y1*2.0;
			double t506 = x0*y3*2.0;
			double t507 = x2*y5*2.0;
			double t502 = t505+t506+t507;
			double t510 = x0*y4*2.0;
			double t511 = x2*y2*2.0;
			double t512 = x1*y5*2.0;
			double t503 = t510+t511+t512;
			double t504 = t501*(1.0/4.0);
			double t508 = t502*t502;
			double t509 = t508*(1.0/4.0);
			double t513 = t503*t503;
			double t514 = t513*(1.0/4.0);
			double t515 = t504+t509+t514;
			double t516 = 1.0/sqrt(t515);
			double t517 = 1.0/pow(t515,3.0/2.0);
			double t518 = 1.0/t515;
			res[0] = -(t516*y4-t503*t517*(t497*y0+t502*y3+t503*y4)*(1.0/4.0))*(t501*t518*(1.0/4.0)-1.0)-t497*t502*t518*(t516*y5-t503*t517*(t497*y3+t502*y1+t503*y5)*(1.0/4.0))*(1.0/4.0)-t497*t503*t518*(t516*y2-t503*t517*(t497*y4+t503*y2+t502*y5)*(1.0/4.0))*(1.0/4.0);
		}
		// gradN2[1]
		{
			double t522 = x1*y1*2.0;
			double t523 = x0*y3*2.0;
			double t524 = x2*y5*2.0;
			double t520 = t522+t523+t524;
			double t527 = x0*y0*2.0;
			double t528 = x1*y3*2.0;
			double t529 = x2*y4*2.0;
			double t521 = t527+t528+t529;
			double t525 = t520*t520;
			double t533 = x0*y4*2.0;
			double t534 = x2*y2*2.0;
			double t535 = x1*y5*2.0;
			double t526 = t533+t534+t535;
			double t530 = t521*t521;
			double t531 = t530*(1.0/4.0);
			double t532 = t525*(1.0/4.0);
			double t536 = t526*t526;
			double t537 = t536*(1.0/4.0);
			double t538 = t531+t532+t537;
			double t539 = 1.0/sqrt(t538);
			double t540 = 1.0/pow(t538,3.0/2.0);
			double t541 = 1.0/t538;
			res[1] = -(t539*y5-t526*t540*(t520*y1+t521*y3+t526*y5)*(1.0/4.0))*(t525*t541*(1.0/4.0)-1.0)-t520*t521*t541*(t539*y4-t526*t540*(t521*y0+t520*y3+t526*y4)*(1.0/4.0))*(1.0/4.0)-t520*t526*t541*(t539*y2-t526*t540*(t520*y5+t521*y4+t526*y2)*(1.0/4.0))*(1.0/4.0);
		}
		// gradN2[2]
		{
			double t546 = x0*y4*2.0;
			double t547 = x2*y2*2.0;
			double t548 = x1*y5*2.0;
			double t543 = t546+t547+t548;
			double t550 = x0*y0*2.0;
			double t551 = x1*y3*2.0;
			double t552 = x2*y4*2.0;
			double t544 = t550+t551+t552;
			double t555 = x1*y1*2.0;
			double t556 = x0*y3*2.0;
			double t557 = x2*y5*2.0;
			double t545 = t555+t556+t557;
			double t549 = t543*t543;
			double t553 = t544*t544;
			double t554 = t553*(1.0/4.0);
			double t558 = t545*t545;
			double t559 = t558*(1.0/4.0);
			double t560 = t549*(1.0/4.0);
			double t561 = t554+t559+t560;
			double t562 = 1.0/sqrt(t561);
			double t563 = 1.0/pow(t561,3.0/2.0);
			double t564 = 1.0/t561;
			res[2] = -(t562*y2-t543*t563*(t543*y2+t544*y4+t545*y5)*(1.0/4.0))*(t549*t564*(1.0/4.0)-1.0)-t543*t544*t564*(t562*y4-t543*t563*(t544*y0+t543*y4+t545*y3)*(1.0/4.0))*(1.0/4.0)-t543*t545*t564*(t562*y5-t543*t563*(t545*y1+t544*y3+t543*y5)*(1.0/4.0))*(1.0/4.0);
		}
		return res ;
	}

};


#endif /* QUADRICGEOMETRYQUANTITIES_H_ */
