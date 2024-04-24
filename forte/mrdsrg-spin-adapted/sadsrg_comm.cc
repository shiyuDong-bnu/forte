/*
 * @BEGIN LICENSE
 *
 * Forte: an open-source plugin to Psi4 (https://github.com/psi4/psi4)
 * that implements a variety of quantum chemistry methods for strongly
 * correlated electrons.
 *
 * Copyright (c) 2012-2024 by its authors (see COPYING, COPYING.LESSER, AUTHORS).
 *
 * The copyrights for code used from other parties are included in
 * the corresponding files.
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see http://www.gnu.org/licenses/.
 *
 * @END LICENSE
 */
#include <algorithm>

#include "psi4/libpsi4util/PsiOutStream.h"

#include "helpers/timer.h"
#include "sadsrg.h"

using namespace psi;

namespace forte {

double SADSRG::H1_T1_C0(BlockedTensor& H1, BlockedTensor& T1, const double& alpha, double& C0) {
    local_timer timer;

    double E = 0.0;
    E += 2.0 * H1["am"] * T1["ma"];

    auto temp = ambit::BlockedTensor::build(tensor_type_, "Temp110", {"aa"});
    temp["uv"] += H1["ev"] * T1["ue"];
    temp["uv"] -= H1["um"] * T1["mv"];

    E += L1_["vu"] * temp["uv"];

    E *= alpha;
    C0 += E;

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H1, T1] -> C0 : %12.3f", timer.get());
    }
    dsrg_time_.add("110", timer.get());
    return E;
}

double SADSRG::H1_T2_C0(BlockedTensor& H1, BlockedTensor& T2, const double& alpha, double& C0) {
    local_timer timer;

    double E = 0.0;
    auto temp = ambit::BlockedTensor::build(tensor_type_, "Temp120", {"aaaa"});
    temp["uvxy"] += H1["ex"] * T2["uvey"];
    temp["uvxy"] -= H1["vm"] * T2["muyx"];

    E += L2_["xyuv"] * temp["uvxy"];

    E *= alpha;
    C0 += E;

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H1, T2] -> C0 : %12.3f", timer.get());
    }
    dsrg_time_.add("120", timer.get());
    return E;
}

double SADSRG::H2_T1_C0(BlockedTensor& H2, BlockedTensor& T1, const double& alpha, double& C0) {
    local_timer timer;

    double E = 0.0;

    auto temp = ambit::BlockedTensor::build(tensor_type_, "Temp120", {"aaaa"});
    temp["uvxy"] += H2["evxy"] * T1["ue"];
    temp["uvxy"] -= H2["uvmy"] * T1["mx"];

    E += L2_["xyuv"] * temp["uvxy"];

    E *= alpha;
    C0 += E;

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T1] -> C0 : %12.3f", timer.get());
    }
    dsrg_time_.add("210", timer.get());
    return E;
}

std::vector<double> SADSRG::H2_T2_C0(BlockedTensor& H2, BlockedTensor& T2, BlockedTensor& S2,
                                     const double& alpha, double& C0) {
    local_timer timer;

    std::vector<double> Eout{0.0, 0.0, 0.0};
    double E = 0.0;

    // [H2, T2] (C_2)^4 from ccvv
    E += H2["efmn"] * S2["mnef"];

    // [H2, T2] (C_2)^4 L1 from cavv
    E += H2["efmu"] * S2["mvef"] * L1_["uv"];

    // [H2, T2] (C_2)^4 L1 from ccav
    E += H2["vemn"] * S2["mnue"] * Eta1_["uv"];

    Eout[0] += E;

    // other terms involving T2 with at least two active indices
    auto Esmall = H2_T2_C0_T2small(H2, T2, S2);

    for (int i = 0; i < 3; ++i) {
        E += Esmall[i];
        Eout[i] += Esmall[i];
        Eout[i] *= alpha;
    }

    // multiply prefactor and copy to C0
    E *= alpha;
    C0 += E;

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T2] -> C0 : %12.3f", timer.get());
    }
    dsrg_time_.add("220", timer.get());
    return Eout;
}

std::vector<double> SADSRG::H2_T2_C0_T2small(BlockedTensor& H2, BlockedTensor& T2,
                                             BlockedTensor& S2) {
    /**
     * Note the following blocks should be available in memory.
     * H2: vvaa, aacc, avca, avac, vaaa, aaca
     * T2: aavv, ccaa, caav, acav, aava, caaa
     * S2: aavv, ccaa, caav, acav, aava, caaa
     */

    double E1 = 0.0, E2 = 0.0, E3 = 0.0;

    // [H2, T2] L1 from aavv
    E1 += 0.25 * H2["efxu"] * S2["yvef"] * L1_["uv"] * L1_["xy"];

    // [H2, T2] L1 from ccaa
    E1 += 0.25 * H2["vymn"] * S2["mnux"] * Eta1_["uv"] * Eta1_["xy"];

    // [H2, T2] L1 from caav
    auto temp = ambit::BlockedTensor::build(tensor_type_, "temp_caav", {"aaaa"});
    temp["uxyv"] += 0.5 * H2["vemx"] * S2["myue"];
    temp["uxyv"] += 0.5 * H2["vexm"] * S2["ymue"];
    E1 += temp["uxyv"] * Eta1_["uv"] * L1_["xy"];

    // [H2, T2] L1 from caaa and aaav
    temp.zero();
    temp.set_name("temp_aaav_caaa");
    temp["uxyv"] += 0.25 * H2["evwx"] * S2["zyeu"] * L1_["wz"];
    temp["uxyv"] += 0.25 * H2["vzmx"] * S2["myuw"] * Eta1_["wz"];
    E1 += temp["uxyv"] * Eta1_["uv"] * L1_["xy"];

    // <[Hbar2, T2]> C_4 (C_2)^2
    temp.zero();
    temp.set_name("temp_H2T2C0_L2");

    // HH
    temp["uvxy"] += 0.5 * H2["uvmn"] * T2["mnxy"];
    temp["uvxy"] += 0.5 * H2["uvmw"] * T2["mzxy"] * L1_["wz"];

    // PP
    temp["uvxy"] += 0.5 * H2["efxy"] * T2["uvef"];
    temp["uvxy"] += 0.5 * H2["ezxy"] * T2["uvew"] * Eta1_["wz"];

    // HP
    temp["uvxy"] += H2["uexm"] * S2["vmye"];
    temp["uvxy"] -= H2["uemx"] * T2["vmye"];
    temp["uvxy"] -= H2["vemx"] * T2["muye"];

    // HP with Gamma1
    temp["uvxy"] += 0.5 * H2["euwx"] * S2["zvey"] * L1_["wz"];
    temp["uvxy"] -= 0.5 * H2["euxw"] * T2["zvey"] * L1_["wz"];
    temp["uvxy"] -= 0.5 * H2["evxw"] * T2["uzey"] * L1_["wz"];

    // HP with Eta1
    temp["uvxy"] += 0.5 * H2["wumx"] * S2["mvzy"] * Eta1_["wz"];
    temp["uvxy"] -= 0.5 * H2["uwmx"] * T2["mvzy"] * Eta1_["wz"];
    temp["uvxy"] -= 0.5 * H2["vwmx"] * T2["muyz"] * Eta1_["wz"];

    E2 += temp["uvxy"] * L2_["uvxy"];

    // <[Hbar2, T2]> C_6 C_2
    if (do_cu3_) {
        if (store_cu3_) {
            timer t("DSRG [H2, T2] L3");
            E3 += H2.block("vaaa")("ewxy") * T2.block("aava")("uvez") * L3_("xyzuwv");
            E3 -= H2.block("aaca")("uvmz") * T2.block("caaa")("mwxy") * L3_("xyzuwv");
        } else {
            // direct algorithm for 3RDM: Alex's trick JCTC 16, 6343–6357 (2020)
            // t_{uvez} v_{ewxy} D_{xyzuwv} = - t_{uvez} v_{ezxy} D_{uvxy}
            //                                + t_{uvez} v_{ewxy} < x^+ y^+ w z^+ v u >

            // - need to transform the integrals to the same orbital basis as active space solver
            // - TODO: maybe we (York) should make the CI vectors consistent at the first place
            ambit::Tensor Tbra, Tket;
            ambit::Tensor Ua = Uactv_.block("aa");

            timer timer_v("DSRG [H2, T2] D3V direct");
            Tbra = H2.block("vaaa").clone();
            Tbra("ewuv") = H2.block("vaaa")("ezxy") * Ua("wz") * Ua("ux") * Ua("vy");
            Tket = T2.block("aava").clone();
            Tket("uvew") = T2.block("aava")("xyez") * Ua("wz") * Ua("ux") * Ua("vy");
            auto E3v_map = as_solver_->compute_complementary_H2caa_overlap(Tbra, Tket);
            timer_v.stop();

            timer timer_c("DSRG [H2, T2] D3C direct");
            Tbra = T2.block("caaa").clone();
            Tbra("mwuv") = T2.block("caaa")("mzxy") * Ua("wz") * Ua("ux") * Ua("vy");
            Tket = H2.block("aaca").clone();
            Tket("uvmw") = H2.block("aaca")("xymz") * Ua("wz") * Ua("ux") * Ua("vy");
            auto E3c_map = as_solver_->compute_complementary_H2caa_overlap(Tbra, Tket);
            timer_c.stop();

            // - 2-RDM contributions
            auto G2 = ambit::BlockedTensor::build(ambit::CoreTensor, "G2", {"aaaa"});
            G2.block("aaaa")("pqrs") = rdms_->SF_G2()("pqrs");

            double E3v = -H2["ezxy"] * T2["uvez"] * G2["xyuv"];
            double E3c = T2["mzxy"] * H2["uvmz"] * G2["xyuv"];

            // - add together
            for (const auto& state_weights : state_to_weights_) {
                const auto& state = state_weights.first;
                const auto& weights = state_weights.second;
                for (size_t i = 0, nroots = weights.size(); i < nroots; ++i) {
                    if (weights[i] < 1.0e-15)
                        continue;
                    E3v += weights[i] * E3v_map[state][i];
                    E3c -= weights[i] * E3c_map[state][i];
                }
            }

            // => spin-free 1- and 2-cumulant contributions <=

            // - virtual contraction
            temp = ambit::BlockedTensor::build(tensor_type_, "temp_va", {"va"});
            temp["ex"] = H2["ewxy"] * L1_["yw"];
            temp["ex"] -= 0.5 * H2["ewyx"] * L1_["yw"];
            E3v -= temp["ex"] * T2["uvez"] * G2["xzuv"];

            temp["eu"] = 0.5 * S2["uvez"] * L1_["zv"];
            E3v -= H2["ewxy"] * temp["eu"] * L2_["xyuw"];

            temp = ambit::BlockedTensor::build(tensor_type_, "temp_vaaa", {"vaaa"});
            temp["ewuy"] = H2["ewxy"] * L1_["xu"];
            E3v -= 0.5 * temp["ewuy"] * S2["uvez"] * L2_["yzwv"];

            temp["ewxu"] = H2["ewxy"] * L1_["yu"];
            E3v += 0.5 * temp["ewxu"] * T2["uvez"] * L2_["xzwv"];
            E3v += 0.5 * temp["ewxv"] * T2["uvez"] * L2_["xzuw"];

            temp["ezxy"] = H2["ewxy"] * L1_["zw"];
            E3v += 0.5 * temp["ezxy"] * T2["uvez"] * G2["xyuv"];

            // - core contraction
            temp = ambit::BlockedTensor::build(tensor_type_, "temp_ac", {"ac"});
            temp["um"] = H2["uvmz"] * L1_["zv"];
            temp["um"] -= 0.5 * H2["vumz"] * L1_["zv"];
            E3c += temp["um"] * T2["mwxy"] * L2_["xyuw"];

            temp["xm"] = S2["mwxy"] * L1_["yw"];
            E3c += 0.5 * H2["uvmz"] * temp["xm"] * G2["xzuv"];

            temp = ambit::BlockedTensor::build(tensor_type_, "temp_caaa", {"caaa"});
            temp["mzxv"] = H2["uvmz"] * L1_["xu"];
            E3c += 0.5 * temp["mzxv"] * S2["mwxy"] * L2_["yzwv"];

            temp["mzuy"] = H2["uvmz"] * L1_["yv"];
            E3c -= 0.5 * temp["mzuy"] * T2["mwxy"] * L2_["xzuw"];
            E3c -= 0.5 * temp["mzux"] * T2["mwxy"] * L2_["yzwu"];

            temp["mwuv"] = H2["uvmz"] * L1_["zw"];
            E3c -= 0.5 * temp["mwuv"] * T2["mwxy"] * G2["xyuv"];

            E3 += E3c + E3v;
        }
    }

    return {E1, E2, E3};
}

void SADSRG::H1_T1_C1(BlockedTensor& H1, BlockedTensor& T1, const double& alpha,
                      BlockedTensor& C1) {
    local_timer timer;

    C1["ip"] += alpha * H1["ap"] * T1["ia"];
    C1["qa"] -= alpha * H1["qi"] * T1["ia"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H1, T1] -> C1 : %12.3f", timer.get());
    }
    dsrg_time_.add("111", timer.get());
}

void SADSRG::H1_T2_C1(BlockedTensor& H1, BlockedTensor& T2, const double& alpha,
                      BlockedTensor& C1) {
    local_timer timer;

    C1["ia"] += 2.0 * alpha * H1["bm"] * T2["imab"];
    C1["ia"] -= alpha * H1["bm"] * T2["miab"];

    C1["ia"] += alpha * H1["bu"] * T2["ivab"] * L1_["uv"];
    C1["ia"] -= 0.5 * alpha * H1["bu"] * T2["viab"] * L1_["uv"];

    C1["ia"] -= alpha * H1["vj"] * T2["ijau"] * L1_["uv"];
    C1["ia"] += 0.5 * alpha * H1["vj"] * T2["jiau"] * L1_["uv"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H1, T2] -> C1 : %12.3f", timer.get());
    }
    dsrg_time_.add("121", timer.get());
}

void SADSRG::H2_T1_C1(BlockedTensor& H2, BlockedTensor& T1, const double& alpha,
                      BlockedTensor& C1) {
    local_timer timer;

    C1["qp"] += 2.0 * alpha * T1["ma"] * H2["qapm"];
    C1["qp"] -= alpha * T1["ma"] * H2["aqpm"];

    C1["qp"] += alpha * T1["xe"] * L1_["yx"] * H2["qepy"];
    C1["qp"] -= 0.5 * alpha * T1["xe"] * L1_["yx"] * H2["eqpy"];

    C1["qp"] -= alpha * T1["mu"] * L1_["uv"] * H2["qvpm"];
    C1["qp"] += 0.5 * alpha * T1["mu"] * L1_["uv"] * H2["vqpm"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T1] -> C1 : %12.3f", timer.get());
    }
    dsrg_time_.add("211", timer.get());
}

void SADSRG::H2_T2_C1(BlockedTensor& H2, BlockedTensor& T2, BlockedTensor& S2, const double& alpha,
                      BlockedTensor& C1) {
    local_timer timer;

    // [Hbar2, T2] (C_2)^3 -> C1 particle contractions
    C1["ir"] += alpha * H2["abrm"] * S2["imab"];

    C1["ir"] += 0.5 * alpha * L1_["uv"] * S2["ivab"] * H2["abru"];

    C1["ir"] += 0.25 * alpha * S2["ijux"] * L1_["xy"] * L1_["uv"] * H2["vyrj"];

    C1["ir"] -= 0.5 * alpha * L1_["uv"] * S2["imub"] * H2["vbrm"];
    C1["ir"] -= 0.5 * alpha * L1_["uv"] * S2["miub"] * H2["bvrm"];

    C1["ir"] -= 0.25 * alpha * S2["iyub"] * L1_["uv"] * L1_["xy"] * H2["vbrx"];
    C1["ir"] -= 0.25 * alpha * S2["iybu"] * L1_["uv"] * L1_["xy"] * H2["bvrx"];

    // [Hbar2, T2] C_4 C_2 2:2 -> C1 ir
    C1["ir"] += 0.5 * alpha * T2["ijxy"] * L2_["xyuv"] * H2["uvrj"];

    C1["ir"] += 0.5 * alpha * H2["aurx"] * S2["ivay"] * L2_["xyuv"];
    C1["ir"] -= 0.5 * alpha * H2["uarx"] * T2["ivay"] * L2_["xyuv"];
    C1["ir"] -= 0.5 * alpha * H2["uarx"] * T2["ivya"] * L2_["xyvu"];

    // [Hbar2, T2] (C_2)^3 -> C1 hole contractions
    C1["pa"] -= alpha * H2["peij"] * S2["ijae"];

    C1["pa"] -= 0.5 * alpha * Eta1_["uv"] * S2["ijau"] * H2["pvij"];

    C1["pa"] -= 0.25 * alpha * S2["vyab"] * Eta1_["uv"] * Eta1_["xy"] * H2["pbux"];

    C1["pa"] += 0.5 * alpha * Eta1_["uv"] * S2["vjae"] * H2["peuj"];
    C1["pa"] += 0.5 * alpha * Eta1_["uv"] * S2["jvae"] * H2["peju"];

    C1["pa"] += 0.25 * alpha * S2["vjax"] * Eta1_["uv"] * Eta1_["xy"] * H2["pyuj"];
    C1["pa"] += 0.25 * alpha * S2["jvax"] * Eta1_["xy"] * Eta1_["uv"] * H2["pyju"];

    // [Hbar2, T2] C_4 C_2 2:2 -> C1 pa
    C1["pa"] -= 0.5 * alpha * L2_["xyuv"] * T2["uvab"] * H2["pbxy"];

    C1["pa"] -= 0.5 * alpha * H2["puix"] * S2["ivay"] * L2_["xyuv"];
    C1["pa"] += 0.5 * alpha * H2["puxi"] * T2["ivay"] * L2_["xyuv"];
    C1["pa"] += 0.5 * alpha * H2["puxi"] * T2["viay"] * L2_["xyvu"];

    // [Hbar2, T2] C_4 C_2 1:3 -> C1
    C1["jb"] += 0.5 * alpha * H2["avxy"] * S2["ujab"] * L2_["xyuv"];

    C1["jb"] -= 0.5 * alpha * H2["uviy"] * S2["ijxb"] * L2_["xyuv"];

    C1["qs"] += alpha * H2["eqxs"] * T2["uvey"] * L2_["xyuv"];
    C1["qs"] -= 0.5 * alpha * H2["eqsx"] * T2["uvey"] * L2_["xyuv"];

    C1["qs"] -= alpha * H2["uqms"] * T2["mvxy"] * L2_["xyuv"];
    C1["qs"] += 0.5 * alpha * H2["uqsm"] * T2["mvxy"] * L2_["xyuv"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T2] -> C1 : %12.3f", timer.get());
    }
    dsrg_time_.add("221", timer.get());
}

void SADSRG::H2_T2_C1_sym(BlockedTensor& H2, BlockedTensor& H2_sym,BlockedTensor& T2, BlockedTensor& S2, const double& alpha,
                  BlockedTensor& C1){
    local_timer timer;
    std::cout<<"Testing function";
    // [Hbar2, T2] (C_2)^3 -> C1 particle contractions
    C1["ir"] += alpha * H2["abrm"] * S2["imab"];

    // C1["ir"] += 0.5 * alpha * L1_["uv"] * S2["ivab"] * H2["abru"];
C1["i,a2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,a4,a5,"] * H2_sym["a4,a5,a2,a9,"];//aaaa
C1["i,v2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,a4,a5,"] * H2_sym["a5,a4,a9,v2,"];//aaav****
C1["i,c2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,a4,a5,"] * H2_sym["a4,a5,c2,a9,"];//aaca
C1["i,a2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,a4,v5,"] * H2_sym["a4,v5,a2,a9,"];//avaa
C1["i,a2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,v4,a5,"] * H2_sym["a5,v4,a9,a2,"];//avaa****
C1["i,c2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,v4,a5,"] * H2_sym["a5,v4,a9,c2,"];//avac****
C1["i,v2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,v4,a5,"] * H2_sym["a5,v4,a9,v2,"];//avav****
C1["i,c2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,a4,v5,"] * H2_sym["a4,v5,c2,a9,"];//avca
C1["i,v2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,a4,v5,"] * H2_sym["a4,v5,v2,a9,"];//avva
C1["i,a2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,v4,v5,"] * H2_sym["v4,v5,a2,a9,"];//vvaa
C1["i,v2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,v4,v5,"] * H2_sym["v5,v4,a9,v2,"];//vvav****
C1["i,c2,"] += 0.5 * alpha * L1_["a9,v,"] * S2["i,v,v4,v5,"] * H2_sym["v4,v5,c2,a9,"];//vvca

    // C1["ir"] += 0.25 * alpha * S2["ijux"] * L1_["xy"] * L1_["uv"] * H2["vyrj"];
C1["i,a2,"] += 0.25 * alpha * S2["i,a1,u,x,"] * L1_["x,a9,"] * L1_["u,a8,"] * H2_sym["a8,a9,a2,a1,"];//aaaa
C1["i,v2,"] += 0.25 * alpha * S2["i,a1,u,x,"] * L1_["x,a9,"] * L1_["u,a8,"] * H2_sym["a9,a8,a1,v2,"];//aaav****
C1["i,c2,"] += 0.25 * alpha * S2["i,a1,u,x,"] * L1_["x,a9,"] * L1_["u,a8,"] * H2_sym["a8,a9,c2,a1,"];//aaca
C1["i,a2,"] += 0.25 * alpha * S2["i,c1,u,x,"] * L1_["x,a9,"] * L1_["u,a8,"] * H2_sym["a9,a8,c1,a2,"];//aaca****
C1["i,c2,"] += 0.25 * alpha * S2["i,c1,u,x,"] * L1_["x,a9,"] * L1_["u,a8,"] * H2_sym["a8,a9,c2,c1,"];//aacc
C1["i,v2,"] += 0.25 * alpha * S2["i,c1,u,x,"] * L1_["x,a9,"] * L1_["u,a8,"] * H2_sym["a9,a8,c1,v2,"];//aacv****

    // C1["ir"] -= 0.5 * alpha * L1_["uv"] * S2["imub"] * H2["vbrm"];
C1["i,a2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["i,c6,u,a5,"] * H2_sym["a5,a8,c6,a2,"];//aaca****
C1["i,c2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["i,c6,u,a5,"] * H2_sym["a8,a5,c2,c6,"];//aacc
C1["i,v2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["i,c6,u,a5,"] * H2_sym["a5,a8,c6,v2,"];//aacv****
C1["i,a2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["i,c6,u,v5,"] * H2_sym["a8,v5,a2,c6,"];//avac
C1["i,c2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["i,c6,u,v5,"] * H2_sym["a8,v5,c2,c6,"];//avcc
C1["i,v2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["i,c6,u,v5,"] * H2_sym["a8,v5,v2,c6,"];//avvc
    // C1["ir"] -= 0.5 * alpha * L1_["uv"] * S2["miub"] * H2["bvrm"];
C1["i,a2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["c6,i,u,a5,"] * H2_sym["a8,a5,c6,a2,"];//aaca****
C1["i,c2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["c6,i,u,a5,"] * H2_sym["a5,a8,c2,c6,"];//aacc
C1["i,v2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["c6,i,u,a5,"] * H2_sym["a8,a5,c6,v2,"];//aacv****
C1["i,a2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["c6,i,u,v5,"] * H2_sym["a8,v5,c6,a2,"];//avca****
C1["i,c2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["c6,i,u,v5,"] * H2_sym["a8,v5,c6,c2,"];//avcc****
C1["i,v2,"] -= 0.5 * alpha * L1_["u,a8,"] * S2["c6,i,u,v5,"] * H2_sym["a8,v5,c6,v2,"];//avcv****

    // C1["ir"] -= 0.25 * alpha * S2["iyub"] * L1_["uv"] * L1_["xy"] * H2["vbrx"];
C1["i,a2,"] -= 0.25 * alpha * S2["i,y,u,a5,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a8,a5,a2,a9,"];//aaaa
C1["i,v2,"] -= 0.25 * alpha * S2["i,y,u,a5,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a5,a8,a9,v2,"];//aaav****
C1["i,c2,"] -= 0.25 * alpha * S2["i,y,u,a5,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a8,a5,c2,a9,"];//aaca
C1["i,a2,"] -= 0.25 * alpha * S2["i,y,u,v5,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a8,v5,a2,a9,"];//avaa
C1["i,c2,"] -= 0.25 * alpha * S2["i,y,u,v5,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a8,v5,c2,a9,"];//avca
C1["i,v2,"] -= 0.25 * alpha * S2["i,y,u,v5,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a8,v5,v2,a9,"];//avva
    // C1["ir"] -= 0.25 * alpha * S2["iybu"] * L1_["uv"] * L1_["xy"] * H2["bvrx"];
C1["i,a2,"] -= 0.25 * alpha * S2["i,y,a5,u,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a5,a8,a2,a9,"];//aaaa
C1["i,v2,"] -= 0.25 * alpha * S2["i,y,a5,u,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a8,a5,a9,v2,"];//aaav****
C1["i,c2,"] -= 0.25 * alpha * S2["i,y,a5,u,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a5,a8,c2,a9,"];//aaca
C1["i,a2,"] -= 0.25 * alpha * S2["i,y,v5,u,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a8,v5,a9,a2,"];//avaa****
C1["i,c2,"] -= 0.25 * alpha * S2["i,y,v5,u,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a8,v5,a9,c2,"];//avac****
C1["i,v2,"] -= 0.25 * alpha * S2["i,y,v5,u,"] * L1_["u,a8,"] * L1_["a9,y,"] * H2_sym["a8,v5,a9,v2,"];//avav****
    // [Hbar2, T2] C_4 C_2 2:2 -> C1 ir
    // C1["ir"] += 0.5 * alpha * T2["ijxy"] * L2_["xyuv"] * H2["uvrj"];
C1["i,a2,"] += 0.5 * alpha * T2["i,a1,x,y,"] * L2_["x,y,a9,a8,"] * H2_sym["a9,a8,a2,a1,"];//aaaa
C1["i,v2,"] += 0.5 * alpha * T2["i,a1,x,y,"] * L2_["x,y,a9,a8,"] * H2_sym["a8,a9,a1,v2,"];//aaav****
C1["i,c2,"] += 0.5 * alpha * T2["i,a1,x,y,"] * L2_["x,y,a9,a8,"] * H2_sym["a9,a8,c2,a1,"];//aaca
C1["i,a2,"] += 0.5 * alpha * T2["i,c1,x,y,"] * L2_["x,y,a9,a8,"] * H2_sym["a8,a9,c1,a2,"];//aaca****
C1["i,c2,"] += 0.5 * alpha * T2["i,c1,x,y,"] * L2_["x,y,a9,a8,"] * H2_sym["a9,a8,c2,c1,"];//aacc
C1["i,v2,"] += 0.5 * alpha * T2["i,c1,x,y,"] * L2_["x,y,a9,a8,"] * H2_sym["a8,a9,c1,v2,"];//aacv****

    // C1["ir"] += 0.5 * alpha * H2["aurx"] * S2["ivay"] * L2_["xyuv"];
C1["i,a2,"] += 0.5 * alpha * H2_sym["a4,a9,a2,a8,"] * S2["i,v,a4,y,"] * L2_["a8,y,a9,v,"];//aaaa
C1["i,v2,"] += 0.5 * alpha * H2_sym["a9,a4,a8,v2,"] * S2["i,v,a4,y,"] * L2_["a8,y,a9,v,"];//aaav****
C1["i,c2,"] += 0.5 * alpha * H2_sym["a4,a9,c2,a8,"] * S2["i,v,a4,y,"] * L2_["a8,y,a9,v,"];//aaca
C1["i,a2,"] += 0.5 * alpha * H2_sym["a9,v4,a8,a2,"] * S2["i,v,v4,y,"] * L2_["a8,y,a9,v,"];//avaa****
C1["i,c2,"] += 0.5 * alpha * H2_sym["a9,v4,a8,c2,"] * S2["i,v,v4,y,"] * L2_["a8,y,a9,v,"];//avac****
C1["i,v2,"] += 0.5 * alpha * H2_sym["a9,v4,a8,v2,"] * S2["i,v,v4,y,"] * L2_["a8,y,a9,v,"];//avav****
    // C1["ir"] -= 0.5 * alpha * H2["uarx"] * T2["ivay"] * L2_["xyuv"];
C1["i,a2,"] -= 0.5 * alpha * H2_sym["a9,a4,a2,a8,"] * T2["i,v,a4,y,"] * L2_["a8,y,a9,v,"];//aaaa
C1["i,v2,"] -= 0.5 * alpha * H2_sym["a4,a9,a8,v2,"] * T2["i,v,a4,y,"] * L2_["a8,y,a9,v,"];//aaav****
C1["i,c2,"] -= 0.5 * alpha * H2_sym["a9,a4,c2,a8,"] * T2["i,v,a4,y,"] * L2_["a8,y,a9,v,"];//aaca
C1["i,a2,"] -= 0.5 * alpha * H2_sym["a9,v4,a2,a8,"] * T2["i,v,v4,y,"] * L2_["a8,y,a9,v,"];//avaa
C1["i,c2,"] -= 0.5 * alpha * H2_sym["a9,v4,c2,a8,"] * T2["i,v,v4,y,"] * L2_["a8,y,a9,v,"];//avca
C1["i,v2,"] -= 0.5 * alpha * H2_sym["a9,v4,v2,a8,"] * T2["i,v,v4,y,"] * L2_["a8,y,a9,v,"];//avva
    // C1["ir"] -= 0.5 * alpha * H2["uarx"] * T2["ivya"] * L2_["xyvu"];
 C1["i,a2,"] -= 0.5 * alpha * H2_sym["a9,a4,a2,a8,"] * T2["i,v,y,a4,"] * L2_["a8,y,v,a9,"];//aaaa
 C1["i,v2,"] -= 0.5 * alpha * H2_sym["a4,a9,a8,v2,"] * T2["i,v,y,a4,"] * L2_["a8,y,v,a9,"];//aaav****
 C1["i,c2,"] -= 0.5 * alpha * H2_sym["a9,a4,c2,a8,"] * T2["i,v,y,a4,"] * L2_["a8,y,v,a9,"];//aaca
 C1["i,a2,"] -= 0.5 * alpha * H2_sym["a9,v4,a2,a8,"] * T2["i,v,y,v4,"] * L2_["a8,y,v,a9,"];//avaa
 C1["i,c2,"] -= 0.5 * alpha * H2_sym["a9,v4,c2,a8,"] * T2["i,v,y,v4,"] * L2_["a8,y,v,a9,"];//avca
 C1["i,v2,"] -= 0.5 * alpha * H2_sym["a9,v4,v2,a8,"] * T2["i,v,y,v4,"] * L2_["a8,y,v,a9,"];//avva

    // [Hbar2, T2] (C_2)^3 -> C1 hole contractions
    // C1["pa"] -= alpha * H2["peij"] * S2["ijae"];
C1["a2,a,"] -= alpha * H2_sym["a2,v7,a0,a1,"] * S2["a0,a1,a,v7,"];//avaa
C1["a2,a,"] -= alpha * H2_sym["a2,v7,a0,c1,"] * S2["a0,c1,a,v7,"];//avac
C1["a2,a,"] -= alpha * H2_sym["a2,v7,c0,a1,"] * S2["c0,a1,a,v7,"];//avca
C1["a2,a,"] -= alpha * H2_sym["a2,v7,c0,c1,"] * S2["c0,c1,a,v7,"];//avcc
C1["c2,a,"] -= alpha * H2_sym["c2,v7,a0,a1,"] * S2["a0,a1,a,v7,"];//cvaa
C1["c2,a,"] -= alpha * H2_sym["c2,v7,a0,c1,"] * S2["a0,c1,a,v7,"];//cvac
C1["c2,a,"] -= alpha * H2_sym["c2,v7,c0,a1,"] * S2["c0,a1,a,v7,"];//cvca
C1["c2,a,"] -= alpha * H2_sym["c2,v7,c0,c1,"] * S2["c0,c1,a,v7,"];//cvcc
C1["v2,a,"] -= alpha * H2_sym["v2,v7,a0,a1,"] * S2["a0,a1,a,v7,"];//vvaa
C1["v2,a,"] -= alpha * H2_sym["v2,v7,c0,a1,"] * S2["c0,a1,a,v7,"];//vvca
C1["v2,a,"] -= alpha * H2_sym["v7,v2,c1,a0,"] * S2["a0,c1,a,v7,"];//vvca****
C1["v2,a,"] -= alpha * H2_sym["v2,v7,c0,c1,"] * S2["c0,c1,a,v7,"];//vvcc

    // C1["pa"] -= 0.5 * alpha * Eta1_["uv"] * S2["ijau"] * H2["pvij"];
C1["a2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["a0,a1,a,u,"] * H2_sym["a2,a8,a0,a1,"];//aaaa
C1["a2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["c0,a1,a,u,"] * H2_sym["a2,a8,c0,a1,"];//aaca
C1["a2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["a0,c1,a,u,"] * H2_sym["a8,a2,c1,a0,"];//aaca****
C1["a2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["c0,c1,a,u,"] * H2_sym["a2,a8,c0,c1,"];//aacc
C1["v2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["a0,a1,a,u,"] * H2_sym["a8,v2,a1,a0,"];//avaa****
C1["v2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["c0,a1,a,u,"] * H2_sym["a8,v2,a1,c0,"];//avac****
C1["v2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["a0,c1,a,u,"] * H2_sym["a8,v2,c1,a0,"];//avca****
C1["v2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["c0,c1,a,u,"] * H2_sym["a8,v2,c1,c0,"];//avcc****
C1["c2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["a0,a1,a,u,"] * H2_sym["c2,a8,a0,a1,"];//caaa
C1["c2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["a0,c1,a,u,"] * H2_sym["c2,a8,a0,c1,"];//caac
C1["c2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["c0,a1,a,u,"] * H2_sym["c2,a8,c0,a1,"];//caca
C1["c2,a,"] -= 0.5 * alpha * Eta1_["u,a8,"] * S2["c0,c1,a,u,"] * H2_sym["c2,a8,c0,c1,"];//cacc

    // C1["pa"] -= 0.25 * alpha * S2["vyab"] * Eta1_["uv"] * Eta1_["xy"] * H2["pbux"];
C1["a2,a,"] -= 0.25 * alpha * S2["v,y,a,a5,"] * Eta1_["a9,v,"] * Eta1_["a8,y,"] * H2_sym["a2,a5,a9,a8,"];//aaaa
C1["a2,a,"] -= 0.25 * alpha * S2["v,y,a,v5,"] * Eta1_["a9,v,"] * Eta1_["a8,y,"] * H2_sym["a2,v5,a9,a8,"];//avaa
C1["v2,a,"] -= 0.25 * alpha * S2["v,y,a,a5,"] * Eta1_["a9,v,"] * Eta1_["a8,y,"] * H2_sym["a5,v2,a8,a9,"];//avaa****
C1["c2,a,"] -= 0.25 * alpha * S2["v,y,a,a5,"] * Eta1_["a9,v,"] * Eta1_["a8,y,"] * H2_sym["c2,a5,a9,a8,"];//caaa
C1["c2,a,"] -= 0.25 * alpha * S2["v,y,a,v5,"] * Eta1_["a9,v,"] * Eta1_["a8,y,"] * H2_sym["c2,v5,a9,a8,"];//cvaa
C1["v2,a,"] -= 0.25 * alpha * S2["v,y,a,v5,"] * Eta1_["a9,v,"] * Eta1_["a8,y,"] * H2_sym["v2,v5,a9,a8,"];//vvaa

    // C1["pa"] += 0.5 * alpha * Eta1_["uv"] * S2["vjae"] * H2["peuj"];
C1["a2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["v,a1,a,v7,"] * H2_sym["a2,v7,a9,a1,"];//avaa
C1["a2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["v,c1,a,v7,"] * H2_sym["a2,v7,a9,c1,"];//avac
C1["c2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["v,a1,a,v7,"] * H2_sym["c2,v7,a9,a1,"];//cvaa
C1["c2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["v,c1,a,v7,"] * H2_sym["c2,v7,a9,c1,"];//cvac
C1["v2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["v,a1,a,v7,"] * H2_sym["v2,v7,a9,a1,"];//vvaa
C1["v2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["v,c1,a,v7,"] * H2_sym["v7,v2,c1,a9,"];//vvca****
    // C1["pa"] += 0.5 * alpha * Eta1_["uv"] * S2["jvae"] * H2["peju"];
C1["a2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["a1,v,a,v7,"] * H2_sym["a2,v7,a1,a9,"];//avaa
C1["a2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["c1,v,a,v7,"] * H2_sym["a2,v7,c1,a9,"];//avca
C1["c2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["a1,v,a,v7,"] * H2_sym["c2,v7,a1,a9,"];//cvaa
C1["c2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["c1,v,a,v7,"] * H2_sym["c2,v7,c1,a9,"];//cvca
C1["v2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["a1,v,a,v7,"] * H2_sym["v2,v7,a1,a9,"];//vvaa
C1["v2,a,"] += 0.5 * alpha * Eta1_["a9,v,"] * S2["c1,v,a,v7,"] * H2_sym["v2,v7,c1,a9,"];//vvca

    // C1["pa"] += 0.25 * alpha * S2["vjax"] * Eta1_["uv"] * Eta1_["xy"] * H2["pyuj"];
C1["a2,a,"] += 0.25 * alpha * S2["v,a9,a,x,"] * Eta1_["a1,v,"] * Eta1_["x,a8,"] * H2_sym["a2,a8,a1,a9,"];//aaaa
C1["a2,a,"] += 0.25 * alpha * S2["v,a9,a,x,"] * Eta1_["c1,v,"] * Eta1_["x,a8,"] * H2_sym["a2,a8,c1,a9,"];//aaca
C1["v2,a,"] += 0.25 * alpha * S2["v,a9,a,x,"] * Eta1_["a1,v,"] * Eta1_["x,a8,"] * H2_sym["a8,v2,a9,a1,"];//avaa****
C1["v2,a,"] += 0.25 * alpha * S2["v,a9,a,x,"] * Eta1_["c1,v,"] * Eta1_["x,a8,"] * H2_sym["a8,v2,a9,c1,"];//avac****
C1["c2,a,"] += 0.25 * alpha * S2["v,a9,a,x,"] * Eta1_["a1,v,"] * Eta1_["x,a8,"] * H2_sym["c2,a8,a1,a9,"];//caaa
C1["c2,a,"] += 0.25 * alpha * S2["v,a9,a,x,"] * Eta1_["c1,v,"] * Eta1_["x,a8,"] * H2_sym["c2,a8,c1,a9,"];//caca
    // C1["pa"] += 0.25 * alpha * S2["jvax"] * Eta1_["xy"] * Eta1_["uv"] * H2["pyju"];
C1["a2,a,"] += 0.25 * alpha * S2["a1,v,a,x,"] * Eta1_["x,a8,"] * Eta1_["a9,v,"] * H2_sym["a2,a8,a1,a9,"];//aaaa
C1["a2,a,"] += 0.25 * alpha * S2["c1,v,a,x,"] * Eta1_["x,a8,"] * Eta1_["a9,v,"] * H2_sym["a2,a8,c1,a9,"];//aaca
C1["v2,a,"] += 0.25 * alpha * S2["a1,v,a,x,"] * Eta1_["x,a8,"] * Eta1_["a9,v,"] * H2_sym["a8,v2,a9,a1,"];//avaa****
C1["v2,a,"] += 0.25 * alpha * S2["c1,v,a,x,"] * Eta1_["x,a8,"] * Eta1_["a9,v,"] * H2_sym["a8,v2,a9,c1,"];//avac****
C1["c2,a,"] += 0.25 * alpha * S2["a1,v,a,x,"] * Eta1_["x,a8,"] * Eta1_["a9,v,"] * H2_sym["c2,a8,a1,a9,"];//caaa
C1["c2,a,"] += 0.25 * alpha * S2["c1,v,a,x,"] * Eta1_["x,a8,"] * Eta1_["a9,v,"] * H2_sym["c2,a8,c1,a9,"];//caca
    // [Hbar2, T2] C_4 C_2 2:2 -> C1 pa
    // C1["pa"] -= 0.5 * alpha * L2_["xyuv"] * T2["uvab"] * H2["pbxy"];
C1["a2,a,"] -= 0.5 * alpha * L2_["a7,a8,u,v,"] * T2["u,v,a,a5,"] * H2_sym["a2,a5,a7,a8,"];//aaaa
C1["a2,a,"] -= 0.5 * alpha * L2_["a7,a8,u,v,"] * T2["u,v,a,v5,"] * H2_sym["a2,v5,a7,a8,"];//avaa
C1["v2,a,"] -= 0.5 * alpha * L2_["a7,a8,u,v,"] * T2["u,v,a,a5,"] * H2_sym["a5,v2,a8,a7,"];//avaa****
C1["c2,a,"] -= 0.5 * alpha * L2_["a7,a8,u,v,"] * T2["u,v,a,a5,"] * H2_sym["c2,a5,a7,a8,"];//caaa
C1["c2,a,"] -= 0.5 * alpha * L2_["a7,a8,u,v,"] * T2["u,v,a,v5,"] * H2_sym["c2,v5,a7,a8,"];//cvaa
C1["v2,a,"] -= 0.5 * alpha * L2_["a7,a8,u,v,"] * T2["u,v,a,v5,"] * H2_sym["v2,v5,a7,a8,"];//vvaa
    // C1["pa"] -= 0.5 * alpha * H2["puix"] * S2["ivay"] * L2_["xyuv"];
C1["a2,a,"] -= 0.5 * alpha * H2_sym["a2,a9,a1,a7,"] * S2["a1,v,a,y,"] * L2_["a7,y,a9,v,"];//aaaa
C1["a2,a,"] -= 0.5 * alpha * H2_sym["a2,a9,c1,a7,"] * S2["c1,v,a,y,"] * L2_["a7,y,a9,v,"];//aaca
C1["v2,a,"] -= 0.5 * alpha * H2_sym["a9,v2,a7,a1,"] * S2["a1,v,a,y,"] * L2_["a7,y,a9,v,"];//avaa****
C1["v2,a,"] -= 0.5 * alpha * H2_sym["a9,v2,a7,c1,"] * S2["c1,v,a,y,"] * L2_["a7,y,a9,v,"];//avac****
C1["c2,a,"] -= 0.5 * alpha * H2_sym["c2,a9,a1,a7,"] * S2["a1,v,a,y,"] * L2_["a7,y,a9,v,"];//caaa
C1["c2,a,"] -= 0.5 * alpha * H2_sym["c2,a9,c1,a7,"] * S2["c1,v,a,y,"] * L2_["a7,y,a9,v,"];//caca
    // C1["pa"] += 0.5 * alpha * H2["puxi"] * T2["ivay"] * L2_["xyuv"];
C1["a2,a,"] += 0.5 * alpha * H2_sym["a2,a9,a7,a1,"] * T2["a1,v,a,y,"] * L2_["a7,y,a9,v,"];//aaaa
C1["a2,a,"] += 0.5 * alpha * H2_sym["a9,a2,c1,a7,"] * T2["c1,v,a,y,"] * L2_["a7,y,a9,v,"];//aaca****
C1["v2,a,"] += 0.5 * alpha * H2_sym["a9,v2,a1,a7,"] * T2["a1,v,a,y,"] * L2_["a7,y,a9,v,"];//avaa****
C1["v2,a,"] += 0.5 * alpha * H2_sym["a9,v2,c1,a7,"] * T2["c1,v,a,y,"] * L2_["a7,y,a9,v,"];//avca****
C1["c2,a,"] += 0.5 * alpha * H2_sym["c2,a9,a7,a1,"] * T2["a1,v,a,y,"] * L2_["a7,y,a9,v,"];//caaa
C1["c2,a,"] += 0.5 * alpha * H2_sym["c2,a9,a7,c1,"] * T2["c1,v,a,y,"] * L2_["a7,y,a9,v,"];//caac
    // C1["pa"] += 0.5 * alpha * H2["puxi"] * T2["viay"] * L2_["xyvu"];
C1["a2,a,"] += 0.5 * alpha * H2_sym["a2,a9,a7,a1,"] * T2["v,a1,a,y,"] * L2_["a7,y,v,a9,"];//aaaa
C1["a2,a,"] += 0.5 * alpha * H2_sym["a9,a2,c1,a7,"] * T2["v,c1,a,y,"] * L2_["a7,y,v,a9,"];//aaca****
C1["v2,a,"] += 0.5 * alpha * H2_sym["a9,v2,a1,a7,"] * T2["v,a1,a,y,"] * L2_["a7,y,v,a9,"];//avaa****
C1["v2,a,"] += 0.5 * alpha * H2_sym["a9,v2,c1,a7,"] * T2["v,c1,a,y,"] * L2_["a7,y,v,a9,"];//avca****
C1["c2,a,"] += 0.5 * alpha * H2_sym["c2,a9,a7,a1,"] * T2["v,a1,a,y,"] * L2_["a7,y,v,a9,"];//caaa
C1["c2,a,"] += 0.5 * alpha * H2_sym["c2,a9,a7,c1,"] * T2["v,c1,a,y,"] * L2_["a7,y,v,a9,"];//caac
    // [Hbar2, T2] C_4 C_2 1:3 -> C1
    // C1["jb"] += 0.5 * alpha * H2["avxy"] * S2["ujab"] * L2_["xyuv"];
C1["j,b,"] += 0.5 * alpha * H2_sym["a4,a9,a7,a8,"] * S2["u,j,a4,b,"] * L2_["a7,a8,u,a9,"];//aaaa
C1["j,b,"] += 0.5 * alpha * H2_sym["a9,v4,a8,a7,"] * S2["u,j,v4,b,"] * L2_["a7,a8,u,a9,"];//avaa****
    // C1["jb"] -= 0.5 * alpha * H2["uviy"] * S2["ijxb"] * L2_["xyuv"];
C1["j,b,"] -= 0.5 * alpha * H2_sym["a7,a9,a1,a8,"] * S2["a1,j,x,b,"] * L2_["x,a8,a7,a9,"];//aaaa
C1["j,b,"] -= 0.5 * alpha * H2_sym["a7,a9,c1,a8,"] * S2["c1,j,x,b,"] * L2_["x,a8,a7,a9,"];//aaca
    // C1["qs"] += alpha * H2["eqxs"] * T2["uvey"] * L2_["xyuv"];
C1["a2,a3,"] += alpha * H2_sym["a2,v9,a3,a8,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//avaa****
C1["a2,c3,"] += alpha * H2_sym["a2,v9,c3,a8,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//avca****
C1["a2,v3,"] += alpha * H2_sym["a2,v9,v3,a8,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//avva****
C1["c2,a3,"] += alpha * H2_sym["c2,v9,a3,a8,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//cvaa****
C1["c2,c3,"] += alpha * H2_sym["c2,v9,c3,a8,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//cvca****
C1["c2,v3,"] += alpha * H2_sym["c2,v9,v3,a8,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//cvva****
C1["v2,a3,"] += alpha * H2_sym["v9,v2,a8,a3,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//vvaa
C1["v2,v3,"] += alpha * H2_sym["v9,v2,a8,v3,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//vvav
C1["v2,c3,"] += alpha * H2_sym["v2,v9,c3,a8,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//vvca****
    // C1["qs"] -= 0.5 * alpha * H2["eqsx"] * T2["uvey"] * L2_["xyuv"];
C1["a2,a3,"] -= 0.5 * alpha * H2_sym["a2,v9,a8,a3,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//avaa****
C1["a2,c3,"] -= 0.5 * alpha * H2_sym["a2,v9,a8,c3,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//avac****
C1["a2,v3,"] -= 0.5 * alpha * H2_sym["a2,v9,a8,v3,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//avav****
C1["c2,a3,"] -= 0.5 * alpha * H2_sym["c2,v9,a8,a3,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//cvaa****
C1["c2,c3,"] -= 0.5 * alpha * H2_sym["c2,v9,a8,c3,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//cvac****
C1["c2,v3,"] -= 0.5 * alpha * H2_sym["c2,v9,a8,v3,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//cvav****
C1["v2,a3,"] -= 0.5 * alpha * H2_sym["v9,v2,a3,a8,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//vvaa
C1["v2,v3,"] -= 0.5 * alpha * H2_sym["v2,v9,a8,v3,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//vvav****
C1["v2,c3,"] -= 0.5 * alpha * H2_sym["v9,v2,c3,a8,"] * T2["u,v,v9,y,"] * L2_["a8,y,u,v,"];//vvca 
    // C1["qs"] -= alpha * H2["uqms"] * T2["mvxy"] * L2_["xyuv"];
 C1["a2,a3,"] -= alpha * H2_sym["a7,a2,c6,a3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//aaca
 C1["a2,c3,"] -= alpha * H2_sym["a7,a2,c6,c3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//aacc
 C1["a2,v3,"] -= alpha * H2_sym["a7,a2,c6,v3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//aacv
 C1["v2,a3,"] -= alpha * H2_sym["a7,v2,c6,a3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//avca
 C1["v2,c3,"] -= alpha * H2_sym["a7,v2,c6,c3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//avcc
 C1["v2,v3,"] -= alpha * H2_sym["a7,v2,c6,v3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//avcv
 C1["c2,a3,"] -= alpha * H2_sym["c2,a7,a3,c6,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//caac****
 C1["c2,c3,"] -= alpha * H2_sym["c2,a7,c3,c6,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//cacc****
 C1["c2,v3,"] -= alpha * H2_sym["c2,a7,v3,c6,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//cavc****
    // C1["qs"] += 0.5 * alpha * H2["uqsm"] * T2["mvxy"] * L2_["xyuv"];
C1["a2,a3,"] += 0.5 * alpha * H2_sym["a2,a7,c6,a3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//aaca****
C1["a2,c3,"] += 0.5 * alpha * H2_sym["a7,a2,c3,c6,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//aacc
C1["a2,v3,"] += 0.5 * alpha * H2_sym["a2,a7,c6,v3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//aacv****
C1["v2,a3,"] += 0.5 * alpha * H2_sym["a7,v2,a3,c6,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//avac
C1["v2,c3,"] += 0.5 * alpha * H2_sym["a7,v2,c3,c6,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//avcc
C1["v2,v3,"] += 0.5 * alpha * H2_sym["a7,v2,v3,c6,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//avvc
C1["c2,a3,"] += 0.5 * alpha * H2_sym["c2,a7,c6,a3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//caca****
C1["c2,c3,"] += 0.5 * alpha * H2_sym["c2,a7,c6,c3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//cacc****
C1["c2,v3,"] += 0.5 * alpha * H2_sym["c2,a7,c6,v3,"] * T2["c6,v,x,y,"] * L2_["x,y,a7,v,"];//cacv****
    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T2] -> C1 : %12.3f", timer.get());
    }
    dsrg_time_.add("221", timer.get());
}

void SADSRG::H1_T2_C2(BlockedTensor& H1, BlockedTensor& T2, const double& alpha,
                      BlockedTensor& C2) {
    local_timer timer;

    C2["ijpb"] += alpha * T2["ijab"] * H1["ap"];
    C2["jibp"] += alpha * T2["ijab"] * H1["ap"];

    C2["qjab"] -= alpha * T2["ijab"] * H1["qi"];
    C2["jqba"] -= alpha * T2["ijab"] * H1["qi"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H1, T2] -> C2 : %12.3f", timer.get());
    }
    dsrg_time_.add("122", timer.get());
}

void SADSRG::H2_T1_C2(BlockedTensor& H2, BlockedTensor& T1, const double& alpha,
                      BlockedTensor& C2) {
    local_timer timer;

    C2["irpq"] += alpha * T1["ia"] * H2["arpq"];
    C2["riqp"] += alpha * T1["ia"] * H2["arpq"];

    C2["rsaq"] -= alpha * T1["ia"] * H2["rsiq"];
    C2["srqa"] -= alpha * T1["ia"] * H2["rsiq"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T1] -> C2 : %12.3f", timer.get());
    }
    dsrg_time_.add("212", timer.get());
}
void SADSRG::H2_T1_C2_sym(BlockedTensor& H2_sym,BlockedTensor& T1, const double& alpha,
                      BlockedTensor& C2) {
    local_timer timer;

    // C2["irpq"] += alpha * T1["ia"] * H2["arpq"];
C2["i,a2,a7,a8,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,a7,a8,"];//aaaa
C2["i,a2,a7,v8,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,a7,v8,"];//aaav
C2["i,a2,v7,a8,"] += alpha * T1["i,a5,"] * H2_sym["a2,a5,a8,v7,"];//aaav****
C2["i,a2,c7,a8,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,c7,a8,"];//aaca
C2["i,a2,a7,c8,"] += alpha * T1["i,a5,"] * H2_sym["a2,a5,c8,a7,"];//aaca****
C2["i,a2,c7,c8,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,c7,c8,"];//aacc
C2["i,a2,c7,v8,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,c7,v8,"];//aacv
C2["i,a2,v7,c8,"] += alpha * T1["i,a5,"] * H2_sym["a2,a5,c8,v7,"];//aacv****
C2["i,a2,v7,v8,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,v7,v8,"];//aavv
C2["i,v2,a7,a8,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,a7,a8,"];//avaa
C2["i,a2,a7,a8,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,a8,a7,"];//avaa****
C2["i,v2,a7,c8,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,a7,c8,"];//avac
C2["i,a2,c7,a8,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,a8,c7,"];//avac****
C2["i,v2,a7,v8,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,a7,v8,"];//avav
C2["i,a2,v7,a8,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,a8,v7,"];//avav****
C2["i,v2,c7,a8,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,c7,a8,"];//avca
C2["i,a2,a7,c8,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,c8,a7,"];//avca****
C2["i,v2,c7,c8,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,c7,c8,"];//avcc
C2["i,a2,c7,c8,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,c8,c7,"];//avcc****
C2["i,v2,c7,v8,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,c7,v8,"];//avcv
C2["i,a2,v7,c8,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,c8,v7,"];//avcv****
C2["i,v2,v7,a8,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,v7,a8,"];//avva
C2["i,a2,a7,v8,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,v8,a7,"];//avva****
C2["i,v2,v7,c8,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,v7,c8,"];//avvc
C2["i,a2,c7,v8,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,v8,c7,"];//avvc****
C2["i,v2,v7,v8,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,v7,v8,"];//avvv
C2["i,a2,v7,v8,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,v8,v7,"];//avvv****
C2["i,c2,a7,a8,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,a8,a7,"];//caaa****
C2["i,c2,c7,a8,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,a8,c7,"];//caac****
C2["i,c2,v7,a8,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,a8,v7,"];//caav****
C2["i,c2,a7,c8,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,c8,a7,"];//caca****
C2["i,c2,c7,c8,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,c8,c7,"];//cacc****
C2["i,c2,v7,c8,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,c8,v7,"];//cacv****
C2["i,c2,a7,v8,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,v8,a7,"];//cava****
C2["i,c2,c7,v8,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,v8,c7,"];//cavc****
C2["i,c2,v7,v8,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,v8,v7,"];//cavv****
C2["i,c2,a7,a8,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,a8,a7,"];//cvaa****
C2["i,c2,c7,a8,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,a8,c7,"];//cvac****
C2["i,c2,v7,a8,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,a8,v7,"];//cvav****
C2["i,c2,a7,c8,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,c8,a7,"];//cvca****
C2["i,c2,c7,c8,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,c8,c7,"];//cvcc****
C2["i,c2,v7,c8,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,c8,v7,"];//cvcv****
C2["i,c2,a7,v8,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,v8,a7,"];//cvva****
C2["i,c2,c7,v8,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,v8,c7,"];//cvvc****
C2["i,c2,v7,v8,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,v8,v7,"];//cvvv****
C2["i,v2,a7,a8,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,a7,a8,"];//vvaa
C2["i,v2,a7,v8,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,a7,v8,"];//vvav
C2["i,v2,v7,a8,"] += alpha * T1["i,v5,"] * H2_sym["v2,v5,a8,v7,"];//vvav****
C2["i,v2,c7,a8,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,c7,a8,"];//vvca
C2["i,v2,a7,c8,"] += alpha * T1["i,v5,"] * H2_sym["v2,v5,c8,a7,"];//vvca****
C2["i,v2,c7,c8,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,c7,c8,"];//vvcc
C2["i,v2,c7,v8,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,c7,v8,"];//vvcv
C2["i,v2,v7,c8,"] += alpha * T1["i,v5,"] * H2_sym["v2,v5,c8,v7,"];//vvcv****
C2["i,v2,v7,v8,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,v7,v8,"];//vvvv
    // C2["riqp"] += alpha * T1["ia"] * H2["arpq"];
C2["a2,i,a8,a7,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,a7,a8,"];//aaaa
C2["a2,i,v8,a7,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,a7,v8,"];//aaav
C2["a2,i,a8,v7,"] += alpha * T1["i,a5,"] * H2_sym["a2,a5,a8,v7,"];//aaav****
C2["a2,i,a8,c7,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,c7,a8,"];//aaca
C2["a2,i,c8,a7,"] += alpha * T1["i,a5,"] * H2_sym["a2,a5,c8,a7,"];//aaca****
C2["a2,i,c8,c7,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,c7,c8,"];//aacc
C2["a2,i,v8,c7,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,c7,v8,"];//aacv
C2["a2,i,c8,v7,"] += alpha * T1["i,a5,"] * H2_sym["a2,a5,c8,v7,"];//aacv****
C2["a2,i,v8,v7,"] += alpha * T1["i,a5,"] * H2_sym["a5,a2,v7,v8,"];//aavv
C2["v2,i,a8,a7,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,a7,a8,"];//avaa
C2["a2,i,a8,a7,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,a8,a7,"];//avaa****
C2["v2,i,c8,a7,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,a7,c8,"];//avac
C2["a2,i,a8,c7,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,a8,c7,"];//avac****
C2["v2,i,v8,a7,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,a7,v8,"];//avav
C2["a2,i,a8,v7,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,a8,v7,"];//avav****
C2["v2,i,a8,c7,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,c7,a8,"];//avca
C2["a2,i,c8,a7,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,c8,a7,"];//avca****
C2["v2,i,c8,c7,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,c7,c8,"];//avcc
C2["a2,i,c8,c7,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,c8,c7,"];//avcc****
C2["v2,i,v8,c7,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,c7,v8,"];//avcv
C2["a2,i,c8,v7,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,c8,v7,"];//avcv****
C2["v2,i,a8,v7,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,v7,a8,"];//avva
C2["a2,i,v8,a7,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,v8,a7,"];//avva****
C2["v2,i,c8,v7,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,v7,c8,"];//avvc
C2["a2,i,v8,c7,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,v8,c7,"];//avvc****
C2["v2,i,v8,v7,"] += alpha * T1["i,a5,"] * H2_sym["a5,v2,v7,v8,"];//avvv
C2["a2,i,v8,v7,"] += alpha * T1["i,v5,"] * H2_sym["a2,v5,v8,v7,"];//avvv****
C2["c2,i,a8,a7,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,a8,a7,"];//caaa****
C2["c2,i,a8,c7,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,a8,c7,"];//caac****
C2["c2,i,a8,v7,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,a8,v7,"];//caav****
C2["c2,i,c8,a7,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,c8,a7,"];//caca****
C2["c2,i,c8,c7,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,c8,c7,"];//cacc****
C2["c2,i,c8,v7,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,c8,v7,"];//cacv****
C2["c2,i,v8,a7,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,v8,a7,"];//cava****
C2["c2,i,v8,c7,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,v8,c7,"];//cavc****
C2["c2,i,v8,v7,"] += alpha * T1["i,a5,"] * H2_sym["c2,a5,v8,v7,"];//cavv****
C2["c2,i,a8,a7,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,a8,a7,"];//cvaa****
C2["c2,i,a8,c7,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,a8,c7,"];//cvac****
C2["c2,i,a8,v7,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,a8,v7,"];//cvav****
C2["c2,i,c8,a7,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,c8,a7,"];//cvca****
C2["c2,i,c8,c7,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,c8,c7,"];//cvcc****
C2["c2,i,c8,v7,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,c8,v7,"];//cvcv****
C2["c2,i,v8,a7,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,v8,a7,"];//cvva****
C2["c2,i,v8,c7,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,v8,c7,"];//cvvc****
C2["c2,i,v8,v7,"] += alpha * T1["i,v5,"] * H2_sym["c2,v5,v8,v7,"];//cvvv****
C2["v2,i,a8,a7,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,a7,a8,"];//vvaa
C2["v2,i,v8,a7,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,a7,v8,"];//vvav
C2["v2,i,a8,v7,"] += alpha * T1["i,v5,"] * H2_sym["v2,v5,a8,v7,"];//vvav****
C2["v2,i,a8,c7,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,c7,a8,"];//vvca
C2["v2,i,c8,a7,"] += alpha * T1["i,v5,"] * H2_sym["v2,v5,c8,a7,"];//vvca****
C2["v2,i,c8,c7,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,c7,c8,"];//vvcc
C2["v2,i,v8,c7,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,c7,v8,"];//vvcv
C2["v2,i,c8,v7,"] += alpha * T1["i,v5,"] * H2_sym["v2,v5,c8,v7,"];//vvcv****
C2["v2,i,v8,v7,"] += alpha * T1["i,v5,"] * H2_sym["v5,v2,v7,v8,"];//vvvv
    // C2["rsaq"] -= alpha * T1["ia"] * H2["rsiq"];
C2["a2,a3,a,a8,"] -= alpha * T1["a9,a,"] * H2_sym["a2,a3,a9,a8,"];//aaaa
C2["a2,a3,a,v8,"] -= alpha * T1["a9,a,"] * H2_sym["a2,a3,a9,v8,"];//aaav
C2["a2,a3,a,a8,"] -= alpha * T1["c9,a,"] * H2_sym["a2,a3,c9,a8,"];//aaca
C2["a2,a3,a,c8,"] -= alpha * T1["a9,a,"] * H2_sym["a3,a2,c8,a9,"];//aaca****
C2["a2,a3,a,c8,"] -= alpha * T1["c9,a,"] * H2_sym["a2,a3,c9,c8,"];//aacc
C2["a2,a3,a,v8,"] -= alpha * T1["c9,a,"] * H2_sym["a2,a3,c9,v8,"];//aacv
C2["a2,v3,a,a8,"] -= alpha * T1["a9,a,"] * H2_sym["a2,v3,a9,a8,"];//avaa
C2["v2,a3,a,a8,"] -= alpha * T1["a9,a,"] * H2_sym["a3,v2,a8,a9,"];//avaa****
C2["a2,v3,a,c8,"] -= alpha * T1["a9,a,"] * H2_sym["a2,v3,a9,c8,"];//avac
C2["v2,a3,a,a8,"] -= alpha * T1["c9,a,"] * H2_sym["a3,v2,a8,c9,"];//avac****
C2["a2,v3,a,v8,"] -= alpha * T1["a9,a,"] * H2_sym["a2,v3,a9,v8,"];//avav
C2["a2,v3,a,a8,"] -= alpha * T1["c9,a,"] * H2_sym["a2,v3,c9,a8,"];//avca
C2["v2,a3,a,c8,"] -= alpha * T1["a9,a,"] * H2_sym["a3,v2,c8,a9,"];//avca****
C2["a2,v3,a,c8,"] -= alpha * T1["c9,a,"] * H2_sym["a2,v3,c9,c8,"];//avcc
C2["v2,a3,a,c8,"] -= alpha * T1["c9,a,"] * H2_sym["a3,v2,c8,c9,"];//avcc****
C2["a2,v3,a,v8,"] -= alpha * T1["c9,a,"] * H2_sym["a2,v3,c9,v8,"];//avcv
C2["v2,a3,a,v8,"] -= alpha * T1["a9,a,"] * H2_sym["a3,v2,v8,a9,"];//avva****
C2["v2,a3,a,v8,"] -= alpha * T1["c9,a,"] * H2_sym["a3,v2,v8,c9,"];//avvc****
C2["c2,a3,a,a8,"] -= alpha * T1["a9,a,"] * H2_sym["c2,a3,a9,a8,"];//caaa
C2["a2,c3,a,a8,"] -= alpha * T1["a9,a,"] * H2_sym["c3,a2,a8,a9,"];//caaa****
C2["c2,a3,a,c8,"] -= alpha * T1["a9,a,"] * H2_sym["c2,a3,a9,c8,"];//caac
C2["a2,c3,a,a8,"] -= alpha * T1["c9,a,"] * H2_sym["c3,a2,a8,c9,"];//caac****
C2["c2,a3,a,v8,"] -= alpha * T1["a9,a,"] * H2_sym["c2,a3,a9,v8,"];//caav
C2["c2,a3,a,a8,"] -= alpha * T1["c9,a,"] * H2_sym["c2,a3,c9,a8,"];//caca
C2["a2,c3,a,c8,"] -= alpha * T1["a9,a,"] * H2_sym["c3,a2,c8,a9,"];//caca****
C2["c2,a3,a,c8,"] -= alpha * T1["c9,a,"] * H2_sym["c2,a3,c9,c8,"];//cacc
C2["a2,c3,a,c8,"] -= alpha * T1["c9,a,"] * H2_sym["c3,a2,c8,c9,"];//cacc****
C2["c2,a3,a,v8,"] -= alpha * T1["c9,a,"] * H2_sym["c2,a3,c9,v8,"];//cacv
C2["a2,c3,a,v8,"] -= alpha * T1["a9,a,"] * H2_sym["c3,a2,v8,a9,"];//cava****
C2["a2,c3,a,v8,"] -= alpha * T1["c9,a,"] * H2_sym["c3,a2,v8,c9,"];//cavc****
C2["c2,c3,a,a8,"] -= alpha * T1["a9,a,"] * H2_sym["c2,c3,a9,a8,"];//ccaa
C2["c2,c3,a,v8,"] -= alpha * T1["a9,a,"] * H2_sym["c2,c3,a9,v8,"];//ccav
C2["c2,c3,a,a8,"] -= alpha * T1["c9,a,"] * H2_sym["c2,c3,c9,a8,"];//ccca
C2["c2,c3,a,c8,"] -= alpha * T1["a9,a,"] * H2_sym["c3,c2,c8,a9,"];//ccca****
C2["c2,c3,a,c8,"] -= alpha * T1["c9,a,"] * H2_sym["c2,c3,c9,c8,"];//cccc
C2["c2,c3,a,v8,"] -= alpha * T1["c9,a,"] * H2_sym["c2,c3,c9,v8,"];//cccv
C2["c2,v3,a,a8,"] -= alpha * T1["a9,a,"] * H2_sym["c2,v3,a9,a8,"];//cvaa
C2["v2,c3,a,a8,"] -= alpha * T1["a9,a,"] * H2_sym["c3,v2,a8,a9,"];//cvaa****
C2["c2,v3,a,c8,"] -= alpha * T1["a9,a,"] * H2_sym["c2,v3,a9,c8,"];//cvac
C2["v2,c3,a,a8,"] -= alpha * T1["c9,a,"] * H2_sym["c3,v2,a8,c9,"];//cvac****
C2["c2,v3,a,v8,"] -= alpha * T1["a9,a,"] * H2_sym["c2,v3,a9,v8,"];//cvav
C2["c2,v3,a,a8,"] -= alpha * T1["c9,a,"] * H2_sym["c2,v3,c9,a8,"];//cvca
C2["v2,c3,a,c8,"] -= alpha * T1["a9,a,"] * H2_sym["c3,v2,c8,a9,"];//cvca****
C2["c2,v3,a,c8,"] -= alpha * T1["c9,a,"] * H2_sym["c2,v3,c9,c8,"];//cvcc
C2["v2,c3,a,c8,"] -= alpha * T1["c9,a,"] * H2_sym["c3,v2,c8,c9,"];//cvcc****
C2["c2,v3,a,v8,"] -= alpha * T1["c9,a,"] * H2_sym["c2,v3,c9,v8,"];//cvcv
C2["v2,c3,a,v8,"] -= alpha * T1["a9,a,"] * H2_sym["c3,v2,v8,a9,"];//cvva****
C2["v2,c3,a,v8,"] -= alpha * T1["c9,a,"] * H2_sym["c3,v2,v8,c9,"];//cvvc****
C2["v2,v3,a,a8,"] -= alpha * T1["a9,a,"] * H2_sym["v2,v3,a9,a8,"];//vvaa
C2["v2,v3,a,v8,"] -= alpha * T1["a9,a,"] * H2_sym["v2,v3,a9,v8,"];//vvav
C2["v2,v3,a,a8,"] -= alpha * T1["c9,a,"] * H2_sym["v2,v3,c9,a8,"];//vvca
C2["v2,v3,a,c8,"] -= alpha * T1["a9,a,"] * H2_sym["v3,v2,c8,a9,"];//vvca****
C2["v2,v3,a,c8,"] -= alpha * T1["c9,a,"] * H2_sym["v2,v3,c9,c8,"];//vvcc
C2["v2,v3,a,v8,"] -= alpha * T1["c9,a,"] * H2_sym["v2,v3,c9,v8,"];//vvcv
    // C2["srqa"] -= alpha * T1["ia"] * H2["rsiq"];
C2["a3,a2,a8,a,"] -= alpha * T1["a9,a,"] * H2_sym["a2,a3,a9,a8,"];//aaaa
C2["a3,a2,v8,a,"] -= alpha * T1["a9,a,"] * H2_sym["a2,a3,a9,v8,"];//aaav
C2["a3,a2,a8,a,"] -= alpha * T1["c9,a,"] * H2_sym["a2,a3,c9,a8,"];//aaca
C2["a3,a2,c8,a,"] -= alpha * T1["a9,a,"] * H2_sym["a3,a2,c8,a9,"];//aaca****
C2["a3,a2,c8,a,"] -= alpha * T1["c9,a,"] * H2_sym["a2,a3,c9,c8,"];//aacc
C2["a3,a2,v8,a,"] -= alpha * T1["c9,a,"] * H2_sym["a2,a3,c9,v8,"];//aacv
C2["v3,a2,a8,a,"] -= alpha * T1["a9,a,"] * H2_sym["a2,v3,a9,a8,"];//avaa
C2["a3,v2,a8,a,"] -= alpha * T1["a9,a,"] * H2_sym["a3,v2,a8,a9,"];//avaa****
C2["v3,a2,c8,a,"] -= alpha * T1["a9,a,"] * H2_sym["a2,v3,a9,c8,"];//avac
C2["a3,v2,a8,a,"] -= alpha * T1["c9,a,"] * H2_sym["a3,v2,a8,c9,"];//avac****
C2["v3,a2,v8,a,"] -= alpha * T1["a9,a,"] * H2_sym["a2,v3,a9,v8,"];//avav
C2["v3,a2,a8,a,"] -= alpha * T1["c9,a,"] * H2_sym["a2,v3,c9,a8,"];//avca
C2["a3,v2,c8,a,"] -= alpha * T1["a9,a,"] * H2_sym["a3,v2,c8,a9,"];//avca****
C2["v3,a2,c8,a,"] -= alpha * T1["c9,a,"] * H2_sym["a2,v3,c9,c8,"];//avcc
C2["a3,v2,c8,a,"] -= alpha * T1["c9,a,"] * H2_sym["a3,v2,c8,c9,"];//avcc****
C2["v3,a2,v8,a,"] -= alpha * T1["c9,a,"] * H2_sym["a2,v3,c9,v8,"];//avcv
C2["a3,v2,v8,a,"] -= alpha * T1["a9,a,"] * H2_sym["a3,v2,v8,a9,"];//avva****
C2["a3,v2,v8,a,"] -= alpha * T1["c9,a,"] * H2_sym["a3,v2,v8,c9,"];//avvc****
C2["a3,c2,a8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c2,a3,a9,a8,"];//caaa
C2["c3,a2,a8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c3,a2,a8,a9,"];//caaa****
C2["a3,c2,c8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c2,a3,a9,c8,"];//caac
C2["c3,a2,a8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c3,a2,a8,c9,"];//caac****
C2["a3,c2,v8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c2,a3,a9,v8,"];//caav
C2["a3,c2,a8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c2,a3,c9,a8,"];//caca
C2["c3,a2,c8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c3,a2,c8,a9,"];//caca****
C2["a3,c2,c8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c2,a3,c9,c8,"];//cacc
C2["c3,a2,c8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c3,a2,c8,c9,"];//cacc****
C2["a3,c2,v8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c2,a3,c9,v8,"];//cacv
C2["c3,a2,v8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c3,a2,v8,a9,"];//cava****
C2["c3,a2,v8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c3,a2,v8,c9,"];//cavc****
C2["c3,c2,a8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c2,c3,a9,a8,"];//ccaa
C2["c3,c2,v8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c2,c3,a9,v8,"];//ccav
C2["c3,c2,a8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c2,c3,c9,a8,"];//ccca
C2["c3,c2,c8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c3,c2,c8,a9,"];//ccca****
C2["c3,c2,c8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c2,c3,c9,c8,"];//cccc
C2["c3,c2,v8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c2,c3,c9,v8,"];//cccv
C2["v3,c2,a8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c2,v3,a9,a8,"];//cvaa
C2["c3,v2,a8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c3,v2,a8,a9,"];//cvaa****
C2["v3,c2,c8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c2,v3,a9,c8,"];//cvac
C2["c3,v2,a8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c3,v2,a8,c9,"];//cvac****
C2["v3,c2,v8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c2,v3,a9,v8,"];//cvav
C2["v3,c2,a8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c2,v3,c9,a8,"];//cvca
C2["c3,v2,c8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c3,v2,c8,a9,"];//cvca****
C2["v3,c2,c8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c2,v3,c9,c8,"];//cvcc
C2["c3,v2,c8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c3,v2,c8,c9,"];//cvcc****
C2["v3,c2,v8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c2,v3,c9,v8,"];//cvcv
C2["c3,v2,v8,a,"] -= alpha * T1["a9,a,"] * H2_sym["c3,v2,v8,a9,"];//cvva****
C2["c3,v2,v8,a,"] -= alpha * T1["c9,a,"] * H2_sym["c3,v2,v8,c9,"];//cvvc****
C2["v3,v2,a8,a,"] -= alpha * T1["a9,a,"] * H2_sym["v2,v3,a9,a8,"];//vvaa
C2["v3,v2,v8,a,"] -= alpha * T1["a9,a,"] * H2_sym["v2,v3,a9,v8,"];//vvav
C2["v3,v2,a8,a,"] -= alpha * T1["c9,a,"] * H2_sym["v2,v3,c9,a8,"];//vvca
C2["v3,v2,c8,a,"] -= alpha * T1["a9,a,"] * H2_sym["v3,v2,c8,a9,"];//vvca****
C2["v3,v2,c8,a,"] -= alpha * T1["c9,a,"] * H2_sym["v2,v3,c9,c8,"];//vvcc
C2["v3,v2,v8,a,"] -= alpha * T1["c9,a,"] * H2_sym["v2,v3,c9,v8,"];//vvcv
    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T1] -> C2 : %12.3f", timer.get());
    }
    dsrg_time_.add("212", timer.get());
}

void SADSRG::H2_T2_C2(BlockedTensor& H2, BlockedTensor& T2, BlockedTensor& S2, const double& alpha,
                      BlockedTensor& C2) {
    local_timer timer;

    // particle-particle contractions
    C2["ijrs"] += alpha * H2["abrs"] * T2["ijab"];

    C2["ijrs"] -= 0.5 * alpha * L1_["xy"] * T2["ijxb"] * H2["ybrs"];
    C2["jisr"] -= 0.5 * alpha * L1_["xy"] * T2["ijxb"] * H2["ybrs"];

    // hole-hole contractions
    C2["pqab"] += alpha * H2["pqij"] * T2["ijab"];

    C2["pqab"] -= 0.5 * alpha * Eta1_["xy"] * T2["yjab"] * H2["pqxj"];
    C2["qpba"] -= 0.5 * alpha * Eta1_["xy"] * T2["yjab"] * H2["pqxj"];

    // hole-particle contractions
    std::vector<std::string> blocks;
    for (const std::string& block : C2.block_labels()) {
        if (block.substr(1, 1) == virt_label_ or block.substr(3, 1) == core_label_)
            continue;
        else
            blocks.push_back(block);
    }

    auto temp = ambit::BlockedTensor::build(tensor_type_, "temp", blocks);
    temp["qjsb"] += alpha * H2["aqms"] * S2["mjab"];
    temp["qjsb"] -= alpha * H2["aqsm"] * T2["mjab"];
    temp["qjsb"] += 0.5 * alpha * L1_["xy"] * S2["yjab"] * H2["aqxs"];
    temp["qjsb"] -= 0.5 * alpha * L1_["xy"] * T2["yjab"] * H2["aqsx"];
    temp["qjsb"] -= 0.5 * alpha * L1_["xy"] * S2["ijxb"] * H2["yqis"];
    temp["qjsb"] += 0.5 * alpha * L1_["xy"] * T2["ijxb"] * H2["yqsi"];

    C2["qjsb"] += temp["qjsb"];
    C2["jqbs"] += temp["qjsb"];

    blocks.clear();
    for (const std::string& block : C2.block_labels()) {
        if (block.substr(0, 1) == virt_label_ or block.substr(3, 1) == core_label_)
            continue;
        else
            blocks.push_back(block);
    }
    // A13
    temp = ambit::BlockedTensor::build(tensor_type_, "temp", blocks);
    temp["jqsb"] -= alpha * H2["aqsm"] * T2["mjba"];
    temp["jqsb"] -= 0.5 * alpha * L1_["xy"] * T2["yjba"] * H2["aqsx"];
    temp["jqsb"] += 0.5 * alpha * L1_["xy"] * T2["ijbx"] * H2["yqsi"];
    // // I want to save the temp  and inputs to compare the value of different terms.
    // temp.save("temp.test",true);
    // L1_.save("L1_.test",true);
    // T2.save("T2.test",true);
    // H2.save("H2.test",true);
    //end
    C2["jqsb"] += temp["jqsb"];
    C2["qjbs"] += temp["jqsb"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T2] -> C2 : %12.3f", timer.get());
    }
    dsrg_time_.add("222", timer.get());
}
void SADSRG::H2_T2_C2_sym(BlockedTensor& H2_sym, BlockedTensor& T2, BlockedTensor& S2, const double& alpha,
                      BlockedTensor& C2) {
    local_timer timer;
   // particle-particle contractions
    // C2["ijrs"] += alpha * H2["abrs"] * T2["ijab"];
C2["i,j,a2,a3,"] += alpha * H2_sym["a5,a4,a2,a3,"] * T2["i,j,a5,a4,"];//aaaa
C2["i,j,a2,v3,"] += alpha * H2_sym["a5,a4,a2,v3,"] * T2["i,j,a5,a4,"];//aaav
C2["i,j,v2,a3,"] += alpha * H2_sym["a4,a5,a3,v2,"] * T2["i,j,a5,a4,"];//aaav****
C2["i,j,c2,a3,"] += alpha * H2_sym["a5,a4,c2,a3,"] * T2["i,j,a5,a4,"];//aaca
C2["i,j,a2,c3,"] += alpha * H2_sym["a4,a5,c3,a2,"] * T2["i,j,a5,a4,"];//aaca****
C2["i,j,c2,c3,"] += alpha * H2_sym["a5,a4,c2,c3,"] * T2["i,j,a5,a4,"];//aacc
C2["i,j,c2,v3,"] += alpha * H2_sym["a5,a4,c2,v3,"] * T2["i,j,a5,a4,"];//aacv
C2["i,j,v2,c3,"] += alpha * H2_sym["a4,a5,c3,v2,"] * T2["i,j,a5,a4,"];//aacv****
C2["i,j,v2,v3,"] += alpha * H2_sym["a5,a4,v2,v3,"] * T2["i,j,a5,a4,"];//aavv
C2["i,j,a2,a3,"] += alpha * H2_sym["a5,v4,a2,a3,"] * T2["i,j,a5,v4,"];//avaa
C2["i,j,a2,a3,"] += alpha * H2_sym["a4,v5,a3,a2,"] * T2["i,j,v5,a4,"];//avaa****
C2["i,j,a2,c3,"] += alpha * H2_sym["a5,v4,a2,c3,"] * T2["i,j,a5,v4,"];//avac
C2["i,j,c2,a3,"] += alpha * H2_sym["a4,v5,a3,c2,"] * T2["i,j,v5,a4,"];//avac****
C2["i,j,a2,v3,"] += alpha * H2_sym["a5,v4,a2,v3,"] * T2["i,j,a5,v4,"];//avav
C2["i,j,v2,a3,"] += alpha * H2_sym["a4,v5,a3,v2,"] * T2["i,j,v5,a4,"];//avav****
C2["i,j,c2,a3,"] += alpha * H2_sym["a5,v4,c2,a3,"] * T2["i,j,a5,v4,"];//avca
C2["i,j,a2,c3,"] += alpha * H2_sym["a4,v5,c3,a2,"] * T2["i,j,v5,a4,"];//avca****
C2["i,j,c2,c3,"] += alpha * H2_sym["a5,v4,c2,c3,"] * T2["i,j,a5,v4,"];//avcc
C2["i,j,c2,c3,"] += alpha * H2_sym["a4,v5,c3,c2,"] * T2["i,j,v5,a4,"];//avcc****
C2["i,j,c2,v3,"] += alpha * H2_sym["a5,v4,c2,v3,"] * T2["i,j,a5,v4,"];//avcv
C2["i,j,v2,c3,"] += alpha * H2_sym["a4,v5,c3,v2,"] * T2["i,j,v5,a4,"];//avcv****
C2["i,j,v2,a3,"] += alpha * H2_sym["a5,v4,v2,a3,"] * T2["i,j,a5,v4,"];//avva
C2["i,j,a2,v3,"] += alpha * H2_sym["a4,v5,v3,a2,"] * T2["i,j,v5,a4,"];//avva****
C2["i,j,v2,c3,"] += alpha * H2_sym["a5,v4,v2,c3,"] * T2["i,j,a5,v4,"];//avvc
C2["i,j,c2,v3,"] += alpha * H2_sym["a4,v5,v3,c2,"] * T2["i,j,v5,a4,"];//avvc****
C2["i,j,v2,v3,"] += alpha * H2_sym["a5,v4,v2,v3,"] * T2["i,j,a5,v4,"];//avvv
C2["i,j,v2,v3,"] += alpha * H2_sym["a4,v5,v3,v2,"] * T2["i,j,v5,a4,"];//avvv****
C2["i,j,a2,a3,"] += alpha * H2_sym["v5,v4,a2,a3,"] * T2["i,j,v5,v4,"];//vvaa
C2["i,j,a2,v3,"] += alpha * H2_sym["v5,v4,a2,v3,"] * T2["i,j,v5,v4,"];//vvav
C2["i,j,v2,a3,"] += alpha * H2_sym["v4,v5,a3,v2,"] * T2["i,j,v5,v4,"];//vvav****
C2["i,j,c2,a3,"] += alpha * H2_sym["v5,v4,c2,a3,"] * T2["i,j,v5,v4,"];//vvca
C2["i,j,a2,c3,"] += alpha * H2_sym["v4,v5,c3,a2,"] * T2["i,j,v5,v4,"];//vvca****
C2["i,j,c2,c3,"] += alpha * H2_sym["v5,v4,c2,c3,"] * T2["i,j,v5,v4,"];//vvcc
C2["i,j,c2,v3,"] += alpha * H2_sym["v5,v4,c2,v3,"] * T2["i,j,v5,v4,"];//vvcv
C2["i,j,v2,c3,"] += alpha * H2_sym["v4,v5,c3,v2,"] * T2["i,j,v5,v4,"];//vvcv****
C2["i,j,v2,v3,"] += alpha * H2_sym["v5,v4,v2,v3,"] * T2["i,j,v5,v4,"];//vvvv
    // C2["ijrs"] -= 0.5 * alpha * L1_["xy"] * T2["ijxb"] * H2["ybrs"];
C2["i,j,a2,a3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,a2,a3,"];//aaaa
C2["i,j,a2,v3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,a2,v3,"];//aaav
C2["i,j,v2,a3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a4,a8,a3,v2,"];//aaav****
C2["i,j,c2,a3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,c2,a3,"];//aaca
C2["i,j,a2,c3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a4,a8,c3,a2,"];//aaca****
C2["i,j,c2,c3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,c2,c3,"];//aacc
C2["i,j,c2,v3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,c2,v3,"];//aacv
C2["i,j,v2,c3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a4,a8,c3,v2,"];//aacv****
C2["i,j,v2,v3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,v2,v3,"];//aavv
C2["i,j,a2,a3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,a2,a3,"];//avaa
C2["i,j,a2,c3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,a2,c3,"];//avac
C2["i,j,a2,v3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,a2,v3,"];//avav
C2["i,j,c2,a3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,c2,a3,"];//avca
C2["i,j,c2,c3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,c2,c3,"];//avcc
C2["i,j,c2,v3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,c2,v3,"];//avcv
C2["i,j,v2,a3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,v2,a3,"];//avva
C2["i,j,v2,c3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,v2,c3,"];//avvc
C2["i,j,v2,v3,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,v2,v3,"];//avvv
    // C2["jisr"] -= 0.5 * alpha * L1_["xy"] * T2["ijxb"] * H2["ybrs"];
C2["j,i,a3,a2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,a2,a3,"];//aaaa
C2["j,i,v3,a2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,a2,v3,"];//aaav
C2["j,i,a3,v2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a4,a8,a3,v2,"];//aaav****
C2["j,i,a3,c2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,c2,a3,"];//aaca
C2["j,i,c3,a2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a4,a8,c3,a2,"];//aaca****
C2["j,i,c3,c2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,c2,c3,"];//aacc
C2["j,i,v3,c2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,c2,v3,"];//aacv
C2["j,i,c3,v2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a4,a8,c3,v2,"];//aacv****
C2["j,i,v3,v2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,a4,"] * H2_sym["a8,a4,v2,v3,"];//aavv
C2["j,i,a3,a2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,a2,a3,"];//avaa
C2["j,i,c3,a2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,a2,c3,"];//avac
C2["j,i,v3,a2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,a2,v3,"];//avav
C2["j,i,a3,c2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,c2,a3,"];//avca
C2["j,i,c3,c2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,c2,c3,"];//avcc
C2["j,i,v3,c2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,c2,v3,"];//avcv
C2["j,i,a3,v2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,v2,a3,"];//avva
C2["j,i,c3,v2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,v2,c3,"];//avvc
C2["j,i,v3,v2,"] -= 0.5 * alpha * L1_["x,a8,"] * T2["i,j,x,v4,"] * H2_sym["a8,v4,v2,v3,"];//avvv
    // hole-hole contractions
    // C2["pqab"] += alpha * H2["pqij"] * T2["ijab"];

C2["a3,a2,a,b,"] += alpha * H2_sym["a3,a2,a9,a1,"] * T2["a9,a1,a,b,"];//aaaa
C2["a3,a2,a,b,"] += alpha * H2_sym["a3,a2,c9,a1,"] * T2["c9,a1,a,b,"];//aaca
C2["a3,a2,a,b,"] += alpha * H2_sym["a2,a3,c1,a9,"] * T2["a9,c1,a,b,"];//aaca****
C2["a3,a2,a,b,"] += alpha * H2_sym["a3,a2,c9,c1,"] * T2["c9,c1,a,b,"];//aacc
C2["a3,v2,a,b,"] += alpha * H2_sym["a3,v2,a9,a1,"] * T2["a9,a1,a,b,"];//avaa
C2["v3,a2,a,b,"] += alpha * H2_sym["a2,v3,a1,a9,"] * T2["a9,a1,a,b,"];//avaa****
C2["a3,v2,a,b,"] += alpha * H2_sym["a3,v2,a9,c1,"] * T2["a9,c1,a,b,"];//avac
C2["v3,a2,a,b,"] += alpha * H2_sym["a2,v3,a1,c9,"] * T2["c9,a1,a,b,"];//avac****
C2["a3,v2,a,b,"] += alpha * H2_sym["a3,v2,c9,a1,"] * T2["c9,a1,a,b,"];//avca
C2["v3,a2,a,b,"] += alpha * H2_sym["a2,v3,c1,a9,"] * T2["a9,c1,a,b,"];//avca****
C2["a3,v2,a,b,"] += alpha * H2_sym["a3,v2,c9,c1,"] * T2["c9,c1,a,b,"];//avcc
C2["v3,a2,a,b,"] += alpha * H2_sym["a2,v3,c1,c9,"] * T2["c9,c1,a,b,"];//avcc****
C2["c3,a2,a,b,"] += alpha * H2_sym["c3,a2,a9,a1,"] * T2["a9,a1,a,b,"];//caaa
C2["a3,c2,a,b,"] += alpha * H2_sym["c2,a3,a1,a9,"] * T2["a9,a1,a,b,"];//caaa****
C2["c3,a2,a,b,"] += alpha * H2_sym["c3,a2,a9,c1,"] * T2["a9,c1,a,b,"];//caac
C2["a3,c2,a,b,"] += alpha * H2_sym["c2,a3,a1,c9,"] * T2["c9,a1,a,b,"];//caac****
C2["c3,a2,a,b,"] += alpha * H2_sym["c3,a2,c9,a1,"] * T2["c9,a1,a,b,"];//caca
C2["a3,c2,a,b,"] += alpha * H2_sym["c2,a3,c1,a9,"] * T2["a9,c1,a,b,"];//caca****
C2["c3,a2,a,b,"] += alpha * H2_sym["c3,a2,c9,c1,"] * T2["c9,c1,a,b,"];//cacc
C2["a3,c2,a,b,"] += alpha * H2_sym["c2,a3,c1,c9,"] * T2["c9,c1,a,b,"];//cacc****
C2["c3,c2,a,b,"] += alpha * H2_sym["c3,c2,a9,a1,"] * T2["a9,a1,a,b,"];//ccaa
C2["c3,c2,a,b,"] += alpha * H2_sym["c3,c2,c9,a1,"] * T2["c9,a1,a,b,"];//ccca
C2["c3,c2,a,b,"] += alpha * H2_sym["c2,c3,c1,a9,"] * T2["a9,c1,a,b,"];//ccca****
C2["c3,c2,a,b,"] += alpha * H2_sym["c3,c2,c9,c1,"] * T2["c9,c1,a,b,"];//cccc
C2["c3,v2,a,b,"] += alpha * H2_sym["c3,v2,a9,a1,"] * T2["a9,a1,a,b,"];//cvaa
C2["v3,c2,a,b,"] += alpha * H2_sym["c2,v3,a1,a9,"] * T2["a9,a1,a,b,"];//cvaa****
C2["c3,v2,a,b,"] += alpha * H2_sym["c3,v2,a9,c1,"] * T2["a9,c1,a,b,"];//cvac
C2["v3,c2,a,b,"] += alpha * H2_sym["c2,v3,a1,c9,"] * T2["c9,a1,a,b,"];//cvac****
C2["c3,v2,a,b,"] += alpha * H2_sym["c3,v2,c9,a1,"] * T2["c9,a1,a,b,"];//cvca
C2["v3,c2,a,b,"] += alpha * H2_sym["c2,v3,c1,a9,"] * T2["a9,c1,a,b,"];//cvca****
C2["c3,v2,a,b,"] += alpha * H2_sym["c3,v2,c9,c1,"] * T2["c9,c1,a,b,"];//cvcc
C2["v3,c2,a,b,"] += alpha * H2_sym["c2,v3,c1,c9,"] * T2["c9,c1,a,b,"];//cvcc****
C2["v3,v2,a,b,"] += alpha * H2_sym["v3,v2,a9,a1,"] * T2["a9,a1,a,b,"];//vvaa
C2["v3,v2,a,b,"] += alpha * H2_sym["v3,v2,c9,a1,"] * T2["c9,a1,a,b,"];//vvca
C2["v3,v2,a,b,"] += alpha * H2_sym["v2,v3,c1,a9,"] * T2["a9,c1,a,b,"];//vvca****
C2["v3,v2,a,b,"] += alpha * H2_sym["v3,v2,c9,c1,"] * T2["c9,c1,a,b,"];//vvcc
    // C2["pqab"] -= 0.5 * alpha * Eta1_["xy"] * T2["yjab"] * H2["pqxj"];
C2["a2,a3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["a2,a3,a7,a1,"];//aaaa
C2["a2,a3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["a3,a2,c1,a7,"];//aaca****
C2["a2,v3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["a2,v3,a7,a1,"];//avaa
C2["v2,a3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["a3,v2,a1,a7,"];//avaa****
C2["a2,v3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["a2,v3,a7,c1,"];//avac
C2["v2,a3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["a3,v2,c1,a7,"];//avca****
C2["c2,a3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c2,a3,a7,a1,"];//caaa
C2["a2,c3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c3,a2,a1,a7,"];//caaa****
C2["c2,a3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c2,a3,a7,c1,"];//caac
C2["a2,c3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c3,a2,c1,a7,"];//caca****
C2["c2,c3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c2,c3,a7,a1,"];//ccaa
C2["c2,c3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c3,c2,c1,a7,"];//ccca****
C2["c2,v3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c2,v3,a7,a1,"];//cvaa
C2["v2,c3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c3,v2,a1,a7,"];//cvaa****
C2["c2,v3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c2,v3,a7,c1,"];//cvac
C2["v2,c3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c3,v2,c1,a7,"];//cvca****
C2["v2,v3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["v2,v3,a7,a1,"];//vvaa
C2["v2,v3,a,b,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["v3,v2,c1,a7,"];//vvca****
    // C2["qpba"] -= 0.5 * alpha * Eta1_["xy"] * T2["yjab"] * H2["pqxj"];
C2["a3,a2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["a2,a3,a7,a1,"];//aaaa
C2["a3,a2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["a3,a2,c1,a7,"];//aaca****
C2["v3,a2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["a2,v3,a7,a1,"];//avaa
C2["a3,v2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["a3,v2,a1,a7,"];//avaa****
C2["v3,a2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["a2,v3,a7,c1,"];//avac
C2["a3,v2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["a3,v2,c1,a7,"];//avca****
C2["a3,c2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c2,a3,a7,a1,"];//caaa
C2["c3,a2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c3,a2,a1,a7,"];//caaa****
C2["a3,c2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c2,a3,a7,c1,"];//caac
C2["c3,a2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c3,a2,c1,a7,"];//caca****
C2["c3,c2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c2,c3,a7,a1,"];//ccaa
C2["c3,c2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c3,c2,c1,a7,"];//ccca****
C2["v3,c2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c2,v3,a7,a1,"];//cvaa
C2["c3,v2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["c3,v2,a1,a7,"];//cvaa****
C2["v3,c2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c2,v3,a7,c1,"];//cvac
C2["c3,v2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["c3,v2,c1,a7,"];//cvca****
C2["v3,v2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,a1,a,b,"] * H2_sym["v2,v3,a7,a1,"];//vvaa
C2["v3,v2,b,a,"] -= 0.5 * alpha * Eta1_["a7,y,"] * T2["y,c1,a,b,"] * H2_sym["v3,v2,c1,a7,"];//vvca****
    // hole-particle contractions
    std::vector<std::string> blocks;
    for (const std::string& block : C2.block_labels()) {
        if (block.substr(1, 1) == virt_label_ or block.substr(3, 1) == core_label_)
            continue;
        else
            blocks.push_back(block);
    }

    auto temp = ambit::BlockedTensor::build(tensor_type_, "temp", blocks);
    // temp["qjsb"] += alpha * H2["aqms"] * S2["mjab"];
temp["a2,j,a3,b,"] += alpha * H2_sym["a5,a2,c6,a3,"] * S2["c6,j,a5,b,"];//aaca
temp["a2,j,c3,b,"] += alpha * H2_sym["a5,a2,c6,c3,"] * S2["c6,j,a5,b,"];//aacc
temp["a2,j,v3,b,"] += alpha * H2_sym["a5,a2,c6,v3,"] * S2["c6,j,a5,b,"];//aacv
temp["a2,j,a3,b,"] += alpha * H2_sym["a2,v5,a3,c6,"] * S2["c6,j,v5,b,"];//avac****
temp["v2,j,a3,b,"] += alpha * H2_sym["a5,v2,c6,a3,"] * S2["c6,j,a5,b,"];//avca
temp["v2,j,c3,b,"] += alpha * H2_sym["a5,v2,c6,c3,"] * S2["c6,j,a5,b,"];//avcc
temp["a2,j,c3,b,"] += alpha * H2_sym["a2,v5,c3,c6,"] * S2["c6,j,v5,b,"];//avcc****
temp["v2,j,v3,b,"] += alpha * H2_sym["a5,v2,c6,v3,"] * S2["c6,j,a5,b,"];//avcv
temp["a2,j,v3,b,"] += alpha * H2_sym["a2,v5,v3,c6,"] * S2["c6,j,v5,b,"];//avvc****
temp["c2,j,a3,b,"] += alpha * H2_sym["c2,a5,a3,c6,"] * S2["c6,j,a5,b,"];//caac****
temp["c2,j,c3,b,"] += alpha * H2_sym["c2,a5,c3,c6,"] * S2["c6,j,a5,b,"];//cacc****
temp["c2,j,v3,b,"] += alpha * H2_sym["c2,a5,v3,c6,"] * S2["c6,j,a5,b,"];//cavc****
temp["c2,j,a3,b,"] += alpha * H2_sym["c2,v5,a3,c6,"] * S2["c6,j,v5,b,"];//cvac****
temp["c2,j,c3,b,"] += alpha * H2_sym["c2,v5,c3,c6,"] * S2["c6,j,v5,b,"];//cvcc****
temp["c2,j,v3,b,"] += alpha * H2_sym["c2,v5,v3,c6,"] * S2["c6,j,v5,b,"];//cvvc****
temp["v2,j,a3,b,"] += alpha * H2_sym["v5,v2,c6,a3,"] * S2["c6,j,v5,b,"];//vvca
temp["v2,j,c3,b,"] += alpha * H2_sym["v5,v2,c6,c3,"] * S2["c6,j,v5,b,"];//vvcc
temp["v2,j,v3,b,"] += alpha * H2_sym["v5,v2,c6,v3,"] * S2["c6,j,v5,b,"];//vvcv
    // temp["qjsb"] -= alpha * H2["aqsm"] * T2["mjab"];
temp["a2,j,a3,b,"] -= alpha * H2_sym["a2,a5,c6,a3,"] * T2["c6,j,a5,b,"];//aaca****
temp["a2,j,c3,b,"] -= alpha * H2_sym["a5,a2,c3,c6,"] * T2["c6,j,a5,b,"];//aacc
temp["a2,j,v3,b,"] -= alpha * H2_sym["a2,a5,c6,v3,"] * T2["c6,j,a5,b,"];//aacv****
temp["v2,j,a3,b,"] -= alpha * H2_sym["a5,v2,a3,c6,"] * T2["c6,j,a5,b,"];//avac
temp["a2,j,a3,b,"] -= alpha * H2_sym["a2,v5,c6,a3,"] * T2["c6,j,v5,b,"];//avca****
temp["v2,j,c3,b,"] -= alpha * H2_sym["a5,v2,c3,c6,"] * T2["c6,j,a5,b,"];//avcc
temp["a2,j,c3,b,"] -= alpha * H2_sym["a2,v5,c6,c3,"] * T2["c6,j,v5,b,"];//avcc****
temp["a2,j,v3,b,"] -= alpha * H2_sym["a2,v5,c6,v3,"] * T2["c6,j,v5,b,"];//avcv****
temp["v2,j,v3,b,"] -= alpha * H2_sym["a5,v2,v3,c6,"] * T2["c6,j,a5,b,"];//avvc
temp["c2,j,a3,b,"] -= alpha * H2_sym["c2,a5,c6,a3,"] * T2["c6,j,a5,b,"];//caca****
temp["c2,j,c3,b,"] -= alpha * H2_sym["c2,a5,c6,c3,"] * T2["c6,j,a5,b,"];//cacc****
temp["c2,j,v3,b,"] -= alpha * H2_sym["c2,a5,c6,v3,"] * T2["c6,j,a5,b,"];//cacv****
temp["c2,j,a3,b,"] -= alpha * H2_sym["c2,v5,c6,a3,"] * T2["c6,j,v5,b,"];//cvca****
temp["c2,j,c3,b,"] -= alpha * H2_sym["c2,v5,c6,c3,"] * T2["c6,j,v5,b,"];//cvcc****
temp["c2,j,v3,b,"] -= alpha * H2_sym["c2,v5,c6,v3,"] * T2["c6,j,v5,b,"];//cvcv****
temp["v2,j,a3,b,"] -= alpha * H2_sym["v2,v5,c6,a3,"] * T2["c6,j,v5,b,"];//vvca****
temp["v2,j,c3,b,"] -= alpha * H2_sym["v5,v2,c3,c6,"] * T2["c6,j,v5,b,"];//vvcc
temp["v2,j,v3,b,"] -= alpha * H2_sym["v2,v5,c6,v3,"] * T2["c6,j,v5,b,"];//vvcv****
    // temp["qjsb"] += 0.5 * alpha * L1_["xy"] * S2["yjab"] * H2["aqxs"];
temp["a3,j,a2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,a5,b,"] * H2_sym["a5,a3,a7,a2,"];//aaaa
temp["a3,j,v2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,a5,b,"] * H2_sym["a5,a3,a7,v2,"];//aaav
temp["a3,j,c2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,a5,b,"] * H2_sym["a3,a5,c2,a7,"];//aaca****
temp["v3,j,a2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,a5,b,"] * H2_sym["a5,v3,a7,a2,"];//avaa
temp["a3,j,a2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,v5,b,"] * H2_sym["a3,v5,a2,a7,"];//avaa****
temp["v3,j,c2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,a5,b,"] * H2_sym["a5,v3,a7,c2,"];//avac
temp["v3,j,v2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,a5,b,"] * H2_sym["a5,v3,a7,v2,"];//avav
temp["a3,j,c2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,v5,b,"] * H2_sym["a3,v5,c2,a7,"];//avca****
temp["a3,j,v2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,v5,b,"] * H2_sym["a3,v5,v2,a7,"];//avva****
temp["c3,j,a2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,a5,b,"] * H2_sym["c3,a5,a2,a7,"];//caaa****
temp["c3,j,c2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,a5,b,"] * H2_sym["c3,a5,c2,a7,"];//caca****
temp["c3,j,v2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,a5,b,"] * H2_sym["c3,a5,v2,a7,"];//cava****
temp["c3,j,a2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,v5,b,"] * H2_sym["c3,v5,a2,a7,"];//cvaa****
temp["c3,j,c2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,v5,b,"] * H2_sym["c3,v5,c2,a7,"];//cvca****
temp["c3,j,v2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,v5,b,"] * H2_sym["c3,v5,v2,a7,"];//cvva****
temp["v3,j,a2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,v5,b,"] * H2_sym["v5,v3,a7,a2,"];//vvaa
temp["v3,j,v2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,v5,b,"] * H2_sym["v5,v3,a7,v2,"];//vvav
temp["v3,j,c2,b,"] += 0.5 * alpha * L1_["a7,y,"] * S2["y,j,v5,b,"] * H2_sym["v3,v5,c2,a7,"];//vvca****
    // temp["qjsb"] -= 0.5 * alpha * L1_["xy"] * T2["yjab"] * H2["aqsx"];
temp["a3,j,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,a5,b,"] * H2_sym["a5,a3,a2,a7,"];//aaaa
temp["a3,j,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,a5,b,"] * H2_sym["a3,a5,a7,v2,"];//aaav****
temp["a3,j,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,a5,b,"] * H2_sym["a5,a3,c2,a7,"];//aaca
temp["v3,j,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,a5,b,"] * H2_sym["a5,v3,a2,a7,"];//avaa
temp["a3,j,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,v5,b,"] * H2_sym["a3,v5,a7,a2,"];//avaa****
temp["a3,j,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,v5,b,"] * H2_sym["a3,v5,a7,c2,"];//avac****
temp["a3,j,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,v5,b,"] * H2_sym["a3,v5,a7,v2,"];//avav****
temp["v3,j,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,a5,b,"] * H2_sym["a5,v3,c2,a7,"];//avca
temp["v3,j,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,a5,b,"] * H2_sym["a5,v3,v2,a7,"];//avva
temp["c3,j,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,a5,b,"] * H2_sym["c3,a5,a7,a2,"];//caaa****
temp["c3,j,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,a5,b,"] * H2_sym["c3,a5,a7,c2,"];//caac****
temp["c3,j,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,a5,b,"] * H2_sym["c3,a5,a7,v2,"];//caav****
temp["c3,j,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,v5,b,"] * H2_sym["c3,v5,a7,a2,"];//cvaa****
temp["c3,j,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,v5,b,"] * H2_sym["c3,v5,a7,c2,"];//cvac****
temp["c3,j,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,v5,b,"] * H2_sym["c3,v5,a7,v2,"];//cvav****
temp["v3,j,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,v5,b,"] * H2_sym["v5,v3,a2,a7,"];//vvaa
temp["v3,j,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,v5,b,"] * H2_sym["v3,v5,a7,v2,"];//vvav****
temp["v3,j,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,v5,b,"] * H2_sym["v5,v3,c2,a7,"];//vvca
    // temp["qjsb"] -= 0.5 * alpha * L1_["xy"] * S2["ijxb"] * H2["yqis"];
temp["a3,j,a2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["a9,j,x,b,"] * H2_sym["a8,a3,a9,a2,"];//aaaa
temp["a3,j,v2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["a9,j,x,b,"] * H2_sym["a8,a3,a9,v2,"];//aaav
temp["a3,j,a2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["c9,j,x,b,"] * H2_sym["a8,a3,c9,a2,"];//aaca
temp["a3,j,c2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["a9,j,x,b,"] * H2_sym["a3,a8,c2,a9,"];//aaca****
temp["a3,j,c2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["c9,j,x,b,"] * H2_sym["a8,a3,c9,c2,"];//aacc
temp["a3,j,v2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["c9,j,x,b,"] * H2_sym["a8,a3,c9,v2,"];//aacv
temp["v3,j,a2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["a9,j,x,b,"] * H2_sym["a8,v3,a9,a2,"];//avaa
temp["v3,j,c2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["a9,j,x,b,"] * H2_sym["a8,v3,a9,c2,"];//avac
temp["v3,j,v2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["a9,j,x,b,"] * H2_sym["a8,v3,a9,v2,"];//avav
temp["v3,j,a2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["c9,j,x,b,"] * H2_sym["a8,v3,c9,a2,"];//avca
temp["v3,j,c2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["c9,j,x,b,"] * H2_sym["a8,v3,c9,c2,"];//avcc
temp["v3,j,v2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["c9,j,x,b,"] * H2_sym["a8,v3,c9,v2,"];//avcv
temp["c3,j,a2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["a9,j,x,b,"] * H2_sym["c3,a8,a2,a9,"];//caaa****
temp["c3,j,a2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["c9,j,x,b,"] * H2_sym["c3,a8,a2,c9,"];//caac****
temp["c3,j,c2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["a9,j,x,b,"] * H2_sym["c3,a8,c2,a9,"];//caca****
temp["c3,j,c2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["c9,j,x,b,"] * H2_sym["c3,a8,c2,c9,"];//cacc****
temp["c3,j,v2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["a9,j,x,b,"] * H2_sym["c3,a8,v2,a9,"];//cava****
temp["c3,j,v2,b,"] -= 0.5 * alpha * L1_["x,a8,"] * S2["c9,j,x,b,"] * H2_sym["c3,a8,v2,c9,"];//cavc****
    // temp["qjsb"] += 0.5 * alpha * L1_["xy"] * T2["ijxb"] * H2["yqsi"];
temp["a3,j,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,x,b,"] * H2_sym["a8,a3,a2,a9,"];//aaaa
temp["a3,j,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,x,b,"] * H2_sym["a3,a8,a9,v2,"];//aaav****
temp["a3,j,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,x,b,"] * H2_sym["a8,a3,c2,a9,"];//aaca
temp["a3,j,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,x,b,"] * H2_sym["a3,a8,c9,a2,"];//aaca****
temp["a3,j,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,x,b,"] * H2_sym["a8,a3,c2,c9,"];//aacc
temp["a3,j,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,x,b,"] * H2_sym["a3,a8,c9,v2,"];//aacv****
temp["v3,j,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,x,b,"] * H2_sym["a8,v3,a2,a9,"];//avaa
temp["v3,j,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,x,b,"] * H2_sym["a8,v3,a2,c9,"];//avac
temp["v3,j,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,x,b,"] * H2_sym["a8,v3,c2,a9,"];//avca
temp["v3,j,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,x,b,"] * H2_sym["a8,v3,c2,c9,"];//avcc
temp["v3,j,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,x,b,"] * H2_sym["a8,v3,v2,a9,"];//avva
temp["v3,j,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,x,b,"] * H2_sym["a8,v3,v2,c9,"];//avvc
temp["c3,j,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,x,b,"] * H2_sym["c3,a8,a9,a2,"];//caaa****
temp["c3,j,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,x,b,"] * H2_sym["c3,a8,a9,c2,"];//caac****
temp["c3,j,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,x,b,"] * H2_sym["c3,a8,a9,v2,"];//caav****
temp["c3,j,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,x,b,"] * H2_sym["c3,a8,c9,a2,"];//caca****
temp["c3,j,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,x,b,"] * H2_sym["c3,a8,c9,c2,"];//cacc****
temp["c3,j,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,x,b,"] * H2_sym["c3,a8,c9,v2,"];//cacv****
    C2["qjsb"] += temp["qjsb"];
    C2["jqbs"] += temp["qjsb"];

    blocks.clear();
    for (const std::string& block : C2.block_labels()) {
        if (block.substr(0, 1) == virt_label_ or block.substr(3, 1) == core_label_)
            continue;
        else
            blocks.push_back(block);
    }
    // A13
    temp = ambit::BlockedTensor::build(tensor_type_, "temp", blocks);
    // temp["jqsb"] -= alpha * H2["aqsm"] * T2["mjba"];
temp["j,a2,a3,b,"] -= alpha * H2_sym["a2,a5,c6,a3,"] * T2["c6,j,b,a5,"];//aaca****
temp["j,a2,c3,b,"] -= alpha * H2_sym["a5,a2,c3,c6,"] * T2["c6,j,b,a5,"];//aacc
temp["j,a2,v3,b,"] -= alpha * H2_sym["a2,a5,c6,v3,"] * T2["c6,j,b,a5,"];//aacv****
temp["j,v2,a3,b,"] -= alpha * H2_sym["a5,v2,a3,c6,"] * T2["c6,j,b,a5,"];//avac
temp["j,a2,a3,b,"] -= alpha * H2_sym["a2,v5,c6,a3,"] * T2["c6,j,b,v5,"];//avca****
temp["j,v2,c3,b,"] -= alpha * H2_sym["a5,v2,c3,c6,"] * T2["c6,j,b,a5,"];//avcc
temp["j,a2,c3,b,"] -= alpha * H2_sym["a2,v5,c6,c3,"] * T2["c6,j,b,v5,"];//avcc****
temp["j,a2,v3,b,"] -= alpha * H2_sym["a2,v5,c6,v3,"] * T2["c6,j,b,v5,"];//avcv****
temp["j,v2,v3,b,"] -= alpha * H2_sym["a5,v2,v3,c6,"] * T2["c6,j,b,a5,"];//avvc
temp["j,c2,a3,b,"] -= alpha * H2_sym["c2,a5,c6,a3,"] * T2["c6,j,b,a5,"];//caca****
temp["j,c2,c3,b,"] -= alpha * H2_sym["c2,a5,c6,c3,"] * T2["c6,j,b,a5,"];//cacc****
temp["j,c2,v3,b,"] -= alpha * H2_sym["c2,a5,c6,v3,"] * T2["c6,j,b,a5,"];//cacv****
temp["j,c2,a3,b,"] -= alpha * H2_sym["c2,v5,c6,a3,"] * T2["c6,j,b,v5,"];//cvca****
temp["j,c2,c3,b,"] -= alpha * H2_sym["c2,v5,c6,c3,"] * T2["c6,j,b,v5,"];//cvcc****
temp["j,c2,v3,b,"] -= alpha * H2_sym["c2,v5,c6,v3,"] * T2["c6,j,b,v5,"];//cvcv****
temp["j,v2,a3,b,"] -= alpha * H2_sym["v2,v5,c6,a3,"] * T2["c6,j,b,v5,"];//vvca****
temp["j,v2,c3,b,"] -= alpha * H2_sym["v5,v2,c3,c6,"] * T2["c6,j,b,v5,"];//vvcc
temp["j,v2,v3,b,"] -= alpha * H2_sym["v2,v5,c6,v3,"] * T2["c6,j,b,v5,"];//vvcv****
    // temp["jqsb"] -= 0.5 * alpha * L1_["xy"] * T2["yjba"] * H2["aqsx"];
temp["j,a3,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,a5,"] * H2_sym["a5,a3,a2,a7,"];//aaaa
temp["j,a3,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,a5,"] * H2_sym["a3,a5,a7,v2,"];//aaav****
temp["j,a3,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,a5,"] * H2_sym["a5,a3,c2,a7,"];//aaca
temp["j,v3,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,a5,"] * H2_sym["a5,v3,a2,a7,"];//avaa
temp["j,a3,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,v5,"] * H2_sym["a3,v5,a7,a2,"];//avaa****
temp["j,a3,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,v5,"] * H2_sym["a3,v5,a7,c2,"];//avac****
temp["j,a3,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,v5,"] * H2_sym["a3,v5,a7,v2,"];//avav****
temp["j,v3,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,a5,"] * H2_sym["a5,v3,c2,a7,"];//avca
temp["j,v3,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,a5,"] * H2_sym["a5,v3,v2,a7,"];//avva
temp["j,c3,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,a5,"] * H2_sym["c3,a5,a7,a2,"];//caaa****
temp["j,c3,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,a5,"] * H2_sym["c3,a5,a7,c2,"];//caac****
temp["j,c3,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,a5,"] * H2_sym["c3,a5,a7,v2,"];//caav****
temp["j,c3,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,v5,"] * H2_sym["c3,v5,a7,a2,"];//cvaa****
temp["j,c3,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,v5,"] * H2_sym["c3,v5,a7,c2,"];//cvac****
temp["j,c3,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,v5,"] * H2_sym["c3,v5,a7,v2,"];//cvav****
temp["j,v3,a2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,v5,"] * H2_sym["v5,v3,a2,a7,"];//vvaa
temp["j,v3,v2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,v5,"] * H2_sym["v3,v5,a7,v2,"];//vvav****
temp["j,v3,c2,b,"] -= 0.5 * alpha * L1_["a7,y,"] * T2["y,j,b,v5,"] * H2_sym["v5,v3,c2,a7,"];//vvca
    // temp["jqsb"] += 0.5 * alpha * L1_["xy"] * T2["ijbx"] * H2["yqsi"];
temp["j,a3,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,b,x,"] * H2_sym["a8,a3,a2,a9,"];//aaaa
temp["j,a3,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,b,x,"] * H2_sym["a3,a8,a9,v2,"];//aaav****
temp["j,a3,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,b,x,"] * H2_sym["a8,a3,c2,a9,"];//aaca
temp["j,a3,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,b,x,"] * H2_sym["a3,a8,c9,a2,"];//aaca****
temp["j,a3,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,b,x,"] * H2_sym["a8,a3,c2,c9,"];//aacc
temp["j,a3,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,b,x,"] * H2_sym["a3,a8,c9,v2,"];//aacv****
temp["j,v3,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,b,x,"] * H2_sym["a8,v3,a2,a9,"];//avaa
temp["j,v3,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,b,x,"] * H2_sym["a8,v3,a2,c9,"];//avac
temp["j,v3,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,b,x,"] * H2_sym["a8,v3,c2,a9,"];//avca
temp["j,v3,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,b,x,"] * H2_sym["a8,v3,c2,c9,"];//avcc
temp["j,v3,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,b,x,"] * H2_sym["a8,v3,v2,a9,"];//avva
temp["j,v3,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,b,x,"] * H2_sym["a8,v3,v2,c9,"];//avvc
temp["j,c3,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,b,x,"] * H2_sym["c3,a8,a9,a2,"];//caaa****
temp["j,c3,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,b,x,"] * H2_sym["c3,a8,a9,c2,"];//caac****
temp["j,c3,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["a9,j,b,x,"] * H2_sym["c3,a8,a9,v2,"];//caav****
temp["j,c3,a2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,b,x,"] * H2_sym["c3,a8,c9,a2,"];//caca****
temp["j,c3,c2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,b,x,"] * H2_sym["c3,a8,c9,c2,"];//cacc****
temp["j,c3,v2,b,"] += 0.5 * alpha * L1_["x,a8,"] * T2["c9,j,b,x,"] * H2_sym["c3,a8,c9,v2,"];//cacv****
    C2["jqsb"] += temp["jqsb"];
    C2["qjbs"] += temp["jqsb"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T2] -> C2 : %12.3f", timer.get());
    }
    dsrg_time_.add("222", timer.get());
}

void SADSRG::V_T1_C0_DF(BlockedTensor& B, BlockedTensor& T1, const double& alpha, double& C0) {
    local_timer timer;

    double E = 0.0;

    auto temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp120", {"Laa"});
    temp["gux"] += B["gex"] * T1["ue"];
    temp["gux"] -= B["gum"] * T1["mx"];

    E += L2_["xyuv"] * temp["gux"] * B["gvy"];

    E *= alpha;
    C0 += E;

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T1] -> C0 : %12.3f", timer.get());
    }
    dsrg_time_.add("210", timer.get());
}

std::vector<double> SADSRG::V_T2_C0_DF(BlockedTensor& B, BlockedTensor& T2, BlockedTensor& S2,
                                       const double& alpha, double& C0) {
    local_timer timer;

    std::vector<double> Eout{0.0, 0.0, 0.0};
    double E = 0.0;

    // [H2, T2] (C_2)^4 from ccvv, cavv, and ccav
    auto temp = ambit::BlockedTensor::build(tensor_type_, "temp_220", {"Lvc"});
    temp["gem"] += B["gfn"] * S2["mnef"];
    temp["gem"] += B["gfu"] * S2["mvef"] * L1_["uv"];
    temp["gem"] += B["gvn"] * S2["nmue"] * Eta1_["uv"];
    E += temp["gem"] * B["gem"];
    Eout[0] += E;

    // form H2 for other blocks that fits memory
    std::vector<std::string> blocks{"aacc", "aaca", "vvaa", "vaaa", "avac", "avca"};
    auto H2 = ambit::BlockedTensor::build(tensor_type_, "temp_H2", blocks);
    H2["abij"] = B["gai"] * B["gbj"];

    auto Esmall = H2_T2_C0_T2small(H2, T2, S2);

    for (int i = 0; i < 3; ++i) {
        Eout[i] += Esmall[i];
    }
    E += Esmall[0] + Esmall[1] + Esmall[2];

    // multiply prefactor and copy to C0
    E *= alpha;
    C0 += E;

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T2] -> C0 : %12.3f", timer.get());
    }
    dsrg_time_.add("220", timer.get());
    return Eout;
}

void SADSRG::V_T1_C1_DF(BlockedTensor& B, BlockedTensor& T1, const double& alpha,
                        BlockedTensor& C1) {
    local_timer timer;

    auto temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp211", {"L"});
    temp["g"] += 2.0 * alpha * T1["ma"] * B["gam"];
    temp["g"] += alpha * T1["xe"] * L1_["yx"] * B["gey"];
    temp["g"] -= alpha * T1["mu"] * L1_["uv"] * B["gvm"];
    C1["qp"] += temp["g"] * B["gqp"];

    temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp211", {"Lgc"});
    temp["gpm"] -= alpha * T1["ma"] * B["gap"];
    temp["gpm"] += 0.5 * alpha * T1["mu"] * L1_["uv"] * B["gvp"];
    C1["qp"] += temp["gpm"] * B["gqm"];

    C1["qp"] -= 0.5 * alpha * T1["xe"] * L1_["yx"] * B["gep"] * B["gqy"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T1] -> C1 : %12.3f", timer.get());
    }
    dsrg_time_.add("211", timer.get());
}

void SADSRG::V_T2_C1_DF(BlockedTensor& B, BlockedTensor& T2, BlockedTensor& S2, const double& alpha,
                        BlockedTensor& C1) {
    local_timer timer;

    // [Hbar2, T2] (C_2)^3 -> C1 particle contractions
    auto temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp221", {"Lhp"});

    temp["gia"] += alpha * B["gbm"] * S2["imab"];

    temp["gia"] += 0.5 * alpha * L1_["uv"] * S2["ivab"] * B["gbu"];

    temp["giv"] += 0.25 * alpha * S2["ijux"] * L1_["xy"] * L1_["uv"] * B["gyj"];

    temp["giv"] -= 0.5 * alpha * L1_["uv"] * S2["imub"] * B["gbm"];
    temp["gia"] -= 0.5 * alpha * L1_["uv"] * S2["miua"] * B["gvm"];

    temp["giv"] -= 0.25 * alpha * S2["iyub"] * L1_["uv"] * L1_["xy"] * B["gbx"];
    temp["gia"] -= 0.25 * alpha * S2["iyau"] * L1_["uv"] * L1_["xy"] * B["gvx"];

    // [Hbar2, T2] C_4 C_2 2:2 -> C1 ir
    temp["giu"] += 0.5 * alpha * T2["ijxy"] * L2_["xyuv"] * B["gvj"];

    temp["gia"] += 0.5 * alpha * B["gux"] * S2["ivay"] * L2_["xyuv"];
    temp["giu"] -= 0.5 * alpha * B["gax"] * T2["ivay"] * L2_["xyuv"];
    temp["giu"] -= 0.5 * alpha * B["gax"] * T2["ivya"] * L2_["xyvu"];

    C1["ir"] += temp["gia"] * B["gar"];

    // [Hbar2, T2] (C_2)^3 -> C1 hole contractions
    temp.zero();

    temp["gia"] -= alpha * B["gej"] * S2["ijae"];

    temp["gia"] -= 0.5 * alpha * Eta1_["uv"] * S2["ijau"] * B["gvj"];

    temp["gua"] -= 0.25 * alpha * S2["vyab"] * Eta1_["uv"] * Eta1_["xy"] * B["gbx"];

    temp["gua"] += 0.5 * alpha * Eta1_["uv"] * S2["vjae"] * B["gej"];
    temp["gia"] += 0.5 * alpha * Eta1_["uv"] * S2["ivae"] * B["geu"];

    temp["gua"] += 0.25 * alpha * S2["vjax"] * Eta1_["uv"] * Eta1_["xy"] * B["gyj"];
    temp["gia"] += 0.25 * alpha * S2["ivax"] * Eta1_["xy"] * Eta1_["uv"] * B["gyu"];

    // [Hbar2, T2] C_4 C_2 2:2 -> C1 pa
    temp["gxa"] -= 0.5 * alpha * L2_["xyuv"] * T2["uvab"] * B["gby"];

    temp["gia"] -= 0.5 * alpha * B["gux"] * S2["ivay"] * L2_["xyuv"];
    temp["gxa"] += 0.5 * alpha * B["gui"] * T2["ivay"] * L2_["xyuv"];
    temp["gxa"] += 0.5 * alpha * B["gui"] * T2["viay"] * L2_["xyvu"];

    C1["pa"] += temp["gia"] * B["gpi"];

    // [Hbar2, T2] C_4 C_2 1:3 -> C1
    temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp221", {"Laa"});
    temp["gxu"] = B["gvy"] * L2_["xyuv"];
    C1["jb"] += 0.5 * alpha * B["gax"] * S2["ujab"] * temp["gxu"];
    C1["jb"] -= 0.5 * alpha * B["gui"] * S2["ijxb"] * temp["gxu"];

    temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp221", {"L"});
    temp["g"] += alpha * B["gex"] * T2["uvey"] * L2_["xyuv"];
    temp["g"] -= alpha * B["gum"] * T2["mvxy"] * L2_["xyuv"];
    C1["qs"] += temp["g"] * B["gqs"];

    C1["qs"] -= 0.5 * alpha * B["ges"] * B["gqx"] * T2["uvey"] * L2_["xyuv"];

    C1["qs"] += 0.5 * alpha * B["gus"] * B["gqm"] * T2["mvxy"] * L2_["xyuv"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T2] -> C1 : %12.3f", timer.get());
    }
    dsrg_time_.add("221", timer.get());
}

void SADSRG::V_T1_C2_DF(BlockedTensor& B, BlockedTensor& T1, const double& alpha,
                        BlockedTensor& C2) {
    local_timer timer;

    C2["irpq"] += alpha * T1["ia"] * B["gap"] * B["grq"];
    C2["riqp"] += alpha * T1["ia"] * B["gap"] * B["grq"];
    C2["rsaq"] -= alpha * T1["ia"] * B["gri"] * B["gsq"];
    C2["srqa"] -= alpha * T1["ia"] * B["gri"] * B["gsq"];

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T1] -> C2 : %12.3f", timer.get());
    }
    dsrg_time_.add("212", timer.get());
}

void SADSRG::V_T2_C2_DF(BlockedTensor& B, BlockedTensor& T2, BlockedTensor& S2, const double& alpha,
                        BlockedTensor& C2) {
    local_timer timer;

    // particle-particle contractions
    C2["ijes"] += batched("e", alpha * B["gae"] * B["gbs"] * T2["ijab"]);
    C2["ijus"] += batched("u", alpha * B["gau"] * B["gbs"] * T2["ijab"]);
    C2["ijms"] += batched("m", alpha * B["gam"] * B["gbs"] * T2["ijab"]);

    std::vector<std::string> C2blocks;
    for (const std::string& block : C2.block_labels()) {
        if (block[0] == virt_label_[0] or block[1] == virt_label_[0])
            continue;
        C2blocks.push_back(block);
    }

    auto temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp222", C2blocks);
    temp["ijes"] += batched("e", L1_["xy"] * T2["ijxb"] * B["gye"] * B["gbs"]);
    temp["ijks"] += L1_["xy"] * T2["ijxb"] * B["gyk"] * B["gbs"];

    C2["ijrs"] -= 0.5 * alpha * temp["ijrs"];
    C2["jisr"] -= 0.5 * alpha * temp["ijrs"];

    // hole-hole contractions
    std::vector<std::string> Vblocks;
    for (const std::string& block : C2.block_labels()) {
        if (block[2] == core_label_[0] or block[3] == core_label_[0])
            continue;
        std::string s = block.substr(0, 2) + "hh";
        if (std::find(Vblocks.begin(), Vblocks.end(), s) == Vblocks.end())
            Vblocks.push_back(s);
    }

    temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp222", Vblocks);
    temp["pqij"] = B["gpi"] * B["gqj"];

    C2["pqab"] += alpha * temp["pqij"] * T2["ijab"];

    C2["pqab"] -= 0.5 * alpha * Eta1_["xy"] * T2["yjab"] * temp["pqxj"];
    C2["qpba"] -= 0.5 * alpha * Eta1_["xy"] * T2["yjab"] * temp["pqxj"];

    // hole-particle contractions
    temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp222", {"Lhp"});
    temp["gjb"] += alpha * B["gam"] * S2["mjab"];
    temp["gjb"] += 0.5 * alpha * L1_["xy"] * S2["yjab"] * B["gax"];
    temp["gjb"] -= 0.5 * alpha * L1_["xy"] * S2["ijxb"] * B["gyi"];

    C2["qjsb"] += temp["gjb"] * B["gqs"];
    C2["jqbs"] += temp["gjb"] * B["gqs"];

    // exchange like terms
    V_T2_C2_DF_PH_X(B, T2, alpha, C2);

    if (print_ > 3) {
        outfile->Printf("\n    Time for [H2, T2] -> C2 : %12.3f", timer.get());
    }
    dsrg_time_.add("222", timer.get());
}

void SADSRG::V_T2_C2_DF_PH_X(BlockedTensor& B, BlockedTensor& T2, const double& alpha,
                             BlockedTensor& C2) {

    std::vector<std::string> qjsb_small, qjsb_large, jqsb_small, jqsb_large;

    for (const std::string& block : C2.block_labels()) {
        auto j = block.substr(1, 1);
        auto b = block.substr(3, 1);

        if (j != virt_label_ and b != core_label_) {
            if (std::count(block.begin(), block.end(), 'v') > 2) {
                qjsb_large.push_back(block);
            } else {
                qjsb_small.push_back(block);
            }
        }

        j = block.substr(0, 1);
        if (j != virt_label_ and b != core_label_) {
            if (std::count(block.begin(), block.end(), 'v') > 2) {
                jqsb_large.push_back(block);
            } else {
                jqsb_small.push_back(block);
            }
        }
    }

    auto temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp222PHX", qjsb_small);
    temp["qjsb"] -= alpha * B["gas"] * B["gqm"] * T2["mjab"];
    temp["qjsb"] -= 0.5 * alpha * L1_["xy"] * T2["yjab"] * B["gas"] * B["gqx"];
    temp["qjsb"] += 0.5 * alpha * L1_["xy"] * T2["ijxb"] * B["gys"] * B["gqi"];

    C2["qjsb"] += temp["qjsb"];
    C2["jqbs"] += temp["qjsb"];

    temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp222PHX", jqsb_small);
    temp["jqsb"] -= alpha * B["gas"] * B["gqm"] * T2["mjba"];
    temp["jqsb"] -= 0.5 * alpha * L1_["xy"] * T2["yjba"] * B["gas"] * B["gqx"];
    temp["jqsb"] += 0.5 * alpha * L1_["xy"] * T2["ijbx"] * B["gys"] * B["gqi"];

    C2["jqsb"] += temp["jqsb"];
    C2["qjbs"] += temp["jqsb"];

    if (!qjsb_large.empty()) {
        C2["e,j,f,v0"] -= batched("e", alpha * B["g,a,f"] * B["g,e,m"] * T2["m,j,a,v0"]);
        C2["j,e,v0,f"] -= batched("e", alpha * B["g,a,f"] * B["g,e,m"] * T2["m,j,a,v0"]);

        temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp222PHX", {"ahpv"});
        temp["xjae"] = L1_["xy"] * T2["yjae"];
        C2["e,j,f,v0"] -= batched("e", 0.5 * alpha * temp["x,j,a,v0"] * B["g,a,f"] * B["g,e,x"]);
        C2["j,e,v0,f"] -= batched("e", 0.5 * alpha * temp["x,j,a,v0"] * B["g,a,f"] * B["g,e,x"]);

        temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp222PHX", {"hhav"});
        temp["ijye"] = L1_["xy"] * T2["ijxe"];
        C2["e,j,f,v0"] += batched("e", 0.5 * alpha * temp["i,j,y,v0"] * B["g,y,f"] * B["g,e,i"]);
        C2["j,e,v0,f"] += batched("e", 0.5 * alpha * temp["i,j,y,v0"] * B["g,y,f"] * B["g,e,i"]);
    }

    if (!jqsb_large.empty()) {
        C2["j,e,f,v0"] -= batched("e", alpha * B["g,a,f"] * B["g,e,m"] * T2["m,j,v0,a"]);
        C2["e,j,v0,f"] -= batched("e", alpha * B["g,a,f"] * B["g,e,m"] * T2["m,j,v0,a"]);

        temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp222PHX", {"ahvp"});
        temp["xjea"] = L1_["xy"] * T2["yjea"];
        C2["j,e,f,v0"] -= batched("e", 0.5 * alpha * temp["x,j,v0,a"] * B["g,a,f"] * B["g,e,x"]);
        C2["e,j,v0,f"] -= batched("e", 0.5 * alpha * temp["x,j,v0,a"] * B["g,a,f"] * B["g,e,x"]);

        temp = ambit::BlockedTensor::build(tensor_type_, "DFtemp222PHX", {"hhva"});
        temp["ijey"] = L1_["xy"] * T2["ijex"];
        C2["j,e,f,v0"] += batched("e", 0.5 * alpha * temp["i,j,v0,y"] * B["g,y,f"] * B["g,e,i"]);
        C2["e,j,v0,f"] += batched("e", 0.5 * alpha * temp["i,j,v0,y"] * B["g,y,f"] * B["g,e,i"]);
    }
}

void SADSRG::H_A_Ca(BlockedTensor& H1, BlockedTensor& H2, BlockedTensor& T1, BlockedTensor& T2,
                    BlockedTensor& S2, const double& alpha, BlockedTensor& C1, BlockedTensor& C2) {
    // set up G2["pqrs"] = 2 * H2["pqrs"] - H2["pqsr"]
    auto G2 = ambit::BlockedTensor::build(tensor_type_, "G2H", {"avac", "aaac", "avaa"});
    G2["pqrs"] = 2.0 * H2["pqrs"] - H2["pqsr"];

    H_A_Ca_small(H1, H2, G2, T1, T2, S2, alpha, C1, C2);

    auto temp = ambit::BlockedTensor::build(ambit::CoreTensor, "tempHACa", {"aa"});
    temp["wz"] += H2["efzm"] * S2["wmef"];
    temp["wz"] -= H2["wemn"] * S2["mnze"];

    C1["uv"] += alpha * temp["uv"];
    C1["vu"] += alpha * temp["uv"];
}

void SADSRG::H_A_Ca_small(BlockedTensor& H1, BlockedTensor& H2, BlockedTensor& G2,
                          BlockedTensor& T1, BlockedTensor& T2, BlockedTensor& S2,
                          const double& alpha, BlockedTensor& C1, BlockedTensor& C2) {
    /**
     * The following blocks should be available in memory:
     * G2: avac, aaac, avaa
     * H2: vvaa, aacc, avca, avac, vaaa, aaca, aaaa
     * T2: aavv, ccaa, caav, acav, aava, caaa, aaaa
     * S2: the same as T2
     */

    auto temp = ambit::BlockedTensor::build(ambit::CoreTensor, "tempHACa", {"aa"});

    H1_T_C1a_smallS(H1, T1, S2, temp);
    H2_T_C1a_smallS(H2, T2, S2, temp);
    H2_T_C1a_smallG(G2, T1, T2, temp);

    C1["uv"] += alpha * temp["uv"];
    C1["vu"] += alpha * temp["uv"];

    temp = ambit::BlockedTensor::build(ambit::CoreTensor, "temp", {"aaaa"});

    H_T_C2a_smallS(H1, H2, T1, T2, S2, temp);

    C2["uvxy"] += alpha * temp["uvxy"];
    C2["xyuv"] += alpha * temp["uvxy"];
}

void SADSRG::H2_T_C1a_smallG(BlockedTensor& G2, BlockedTensor& T1, BlockedTensor& T2,
                             BlockedTensor& C1) {
    /**
     * The following blocks should be available in memory:
     * G2: avac, aaac, avaa
     * T2: aava, caaa
     */

    C1["uv"] += T1["ma"] * G2["uavm"];
    C1["uv"] += 0.5 * T1["xe"] * L1_["yx"] * G2["uevy"];
    C1["uv"] -= 0.5 * T1["mx"] * L1_["xy"] * G2["uyvm"];

    C1["wz"] += 0.5 * G2["wezx"] * T2["uvey"] * L2_["xyuv"];
    C1["wz"] -= 0.5 * G2["wuzm"] * T2["mvxy"] * L2_["xyuv"];
}

void SADSRG::H1_T_C1a_smallS(BlockedTensor& H1, BlockedTensor& T1, BlockedTensor& S2,
                             BlockedTensor& C1) {
    // The following blocks should be available in memory:
    // S2: aavv, ccaa, caav, acav, aava, caaa, aaaa

    C1["uv"] += H1["ev"] * T1["ue"];
    C1["uv"] -= H1["um"] * T1["mv"];

    // C1["uv"] += H1["bm"] * S2["umvb"];
    C1["uv"] += H1["em"] * S2["umve"];
    C1["uv"] += H1["xm"] * S2["muxv"];

    // C1["uv"] += 0.5 * H1["bx"] * S2["uyvb"] * L1_["xy"];
    C1["uv"] += 0.5 * H1["ex"] * S2["yuev"] * L1_["xy"];
    C1["uv"] += 0.5 * H1["wx"] * S2["uyvw"] * L1_["xy"];

    // C1["uv"] -= 0.5 * H1["yj"] * S2["ujvx"] * L1_["xy"];
    C1["uv"] -= 0.5 * H1["ym"] * S2["muxv"] * L1_["xy"];
    C1["uv"] -= 0.5 * H1["yw"] * S2["uwvx"] * L1_["xy"];
}

void SADSRG::H2_T_C1a_smallS(BlockedTensor& H2, BlockedTensor& T2, BlockedTensor& S2,
                             BlockedTensor& C1) {
    /**
     * The following blocks should be available in memory:
     * H2: vvaa, aacc, avca, avac, vaaa, aaca, aaaa
     * T2: aavv, ccaa, caav, acav, aava, caaa, aaaa
     * S2: the same as T2
     */

    C1["wz"] += H2["uemz"] * S2["mwue"];
    C1["wz"] += H2["uezm"] * S2["wmue"];
    C1["wz"] += H2["vumz"] * S2["mwvu"];

    C1["wz"] -= H2["wemu"] * S2["muze"];
    C1["wz"] -= H2["weum"] * S2["umze"];
    C1["wz"] -= H2["ewvu"] * S2["vuez"];

    auto temp = ambit::BlockedTensor::build(ambit::CoreTensor, "temp", {"aaaa"});

    // temp["wzuv"] += 0.5 * S2["wvab"] * H2["abzu"];
    temp["wzuv"] += 0.5 * S2["wvef"] * H2["efzu"];
    temp["wzuv"] += 0.5 * S2["wvex"] * H2["exzu"];
    temp["wzuv"] += 0.5 * S2["vwex"] * H2["exuz"];
    temp["wzuv"] += 0.5 * S2["wvxy"] * H2["xyzu"];

    // temp["wzuv"] -= 0.5 * S2["wmub"] * H2["vbzm"];
    temp["wzuv"] -= 0.5 * S2["wmue"] * H2["vezm"];
    temp["wzuv"] -= 0.5 * S2["mwxu"] * H2["xvmz"];

    // temp["wzuv"] -= 0.5 * S2["mwub"] * H2["bvzm"];
    temp["wzuv"] -= 0.5 * S2["mwue"] * H2["vemz"];
    temp["wzuv"] -= 0.5 * S2["mwux"] * H2["vxmz"];

    temp["wzuv"] += 0.25 * S2["jwxu"] * L1_["xy"] * H2["yvjz"];
    temp["wzuv"] -= 0.25 * S2["ywbu"] * L1_["xy"] * H2["bvxz"];
    temp["wzuv"] -= 0.25 * S2["wybu"] * L1_["xy"] * H2["bvzx"];
    C1["wz"] += temp["wzuv"] * L1_["uv"];

    temp.zero();

    // temp["wzuv"] -= 0.5 * S2["ijzu"] * H2["wvij"];
    temp["wzuv"] -= 0.5 * S2["mnzu"] * H2["wvmn"];
    temp["wzuv"] -= 0.5 * S2["mxzu"] * H2["wvmx"];
    temp["wzuv"] -= 0.5 * S2["mxuz"] * H2["vwmx"];
    temp["wzuv"] -= 0.5 * S2["xyzu"] * H2["wvxy"];

    // temp["wzuv"] += 0.5 * S2["vjze"] * H2["weuj"];
    temp["wzuv"] += 0.5 * S2["vmze"] * H2["weum"];
    temp["wzuv"] += 0.5 * S2["xvez"] * H2["ewxu"];

    // temp["wzuv"] += 0.5 * S2["jvze"] * H2["weju"];
    temp["wzuv"] += 0.5 * S2["mvze"] * H2["wemu"];
    temp["wzuv"] += 0.5 * S2["vxez"] * H2["ewux"];

    temp["wzuv"] -= 0.25 * S2["yvbz"] * Eta1_["xy"] * H2["bwxu"];
    temp["wzuv"] += 0.25 * S2["jvxz"] * Eta1_["xy"] * H2["ywju"];
    temp["wzuv"] += 0.25 * S2["jvzx"] * Eta1_["xy"] * H2["wyju"];
    C1["wz"] += temp["wzuv"] * Eta1_["uv"];

    C1["wz"] += 0.5 * H2["vujz"] * T2["jwyx"] * L2_["xyuv"];
    C1["wz"] += 0.5 * H2["auzx"] * S2["wvay"] * L2_["xyuv"];
    C1["wz"] -= 0.5 * H2["auxz"] * T2["wvay"] * L2_["xyuv"];
    C1["wz"] -= 0.5 * H2["auxz"] * T2["vway"] * L2_["xyvu"];

    C1["wz"] -= 0.5 * H2["bwyx"] * T2["vubz"] * L2_["xyuv"];
    C1["wz"] -= 0.5 * H2["wuix"] * S2["ivzy"] * L2_["xyuv"];
    C1["wz"] += 0.5 * H2["uwix"] * T2["ivzy"] * L2_["xyuv"];
    C1["wz"] += 0.5 * H2["uwix"] * T2["ivyz"] * L2_["xyvu"];

    C1["wz"] += 0.5 * H2["avxy"] * S2["uwaz"] * L2_["xyuv"];
    C1["wz"] -= 0.5 * H2["uviy"] * S2["iwxz"] * L2_["xyuv"];
}

void SADSRG::H_T_C2a_smallS(BlockedTensor& H1, BlockedTensor& H2, BlockedTensor& T1,
                            BlockedTensor& T2, BlockedTensor& S2, BlockedTensor& C2) {
    /**
     * The following blocks should be available in memory:
     * H2: vvaa, aacc, avca, avac, vaaa, aaca, aaaa
     * T2: aavv, ccaa, caav, acav, aava, caaa, aaaa
     * S2: the same as T2
     */

    // C2["uvxy"] += H2["abxy"] * T2["uvab"];
    C2["uvxy"] += H2["efxy"] * T2["uvef"];
    C2["uvxy"] += H2["wzxy"] * T2["uvwz"];
    C2["uvxy"] += H2["ewxy"] * T2["uvew"];
    C2["uvxy"] += H2["ewyx"] * T2["vuew"];

    // C2["uvxy"] += H2["uvij"] * T2["ijxy"];
    C2["uvxy"] += H2["uvmn"] * T2["mnxy"];
    C2["uvxy"] += H2["uvwz"] * T2["wzxy"];
    C2["uvxy"] += H2["vumw"] * T2["mwyx"];
    C2["uvxy"] += H2["uvmw"] * T2["mwxy"];

    auto temp = ambit::BlockedTensor::build(ambit::CoreTensor, "temp", {"aaaa"});
    temp["uvxy"] += H1["ax"] * T2["uvay"];
    temp["uvxy"] -= H1["ui"] * T2["ivxy"];
    temp["uvxy"] += T1["ua"] * H2["avxy"];
    temp["uvxy"] -= T1["ix"] * H2["uviy"];

    temp["uvxy"] -= 0.5 * L1_["wz"] * T2["vuaw"] * H2["azyx"];
    temp["uvxy"] -= 0.5 * Eta1_["wz"] * T2["izyx"] * H2["vuiw"];

    // temp["uvxy"] += H2["aumx"] * S2["mvay"];
    temp["uvxy"] += H2["uexm"] * S2["vmye"];
    temp["uvxy"] += H2["wumx"] * S2["mvwy"];

    temp["uvxy"] += 0.5 * L1_["wz"] * S2["zvay"] * H2["auwx"];
    temp["uvxy"] -= 0.5 * L1_["wz"] * S2["ivwy"] * H2["zuix"];

    // temp["uvxy"] -= H2["auxm"] * T2["mvay"];
    temp["uvxy"] -= H2["uemx"] * T2["vmye"];
    temp["uvxy"] -= H2["uwmx"] * T2["mvwy"];

    temp["uvxy"] -= 0.5 * L1_["wz"] * T2["zvay"] * H2["auxw"];
    temp["uvxy"] += 0.5 * L1_["wz"] * T2["ivwy"] * H2["uzix"];

    // temp["uvxy"] -= H2["avxm"] * T2["muya"];
    temp["uvxy"] -= H2["vemx"] * T2["muye"];
    temp["uvxy"] -= H2["vwmx"] * T2["muyw"];

    temp["uvxy"] -= 0.5 * L1_["wz"] * T2["uzay"] * H2["avxw"];
    temp["uvxy"] += 0.5 * L1_["wz"] * T2["iuyw"] * H2["vzix"];

    C2["uvxy"] += temp["uvxy"];
    C2["vuyx"] += temp["uvxy"];
}

void SADSRG::H1d_A1_C1ph(BlockedTensor& H1, BlockedTensor& T1, const double& alpha,
                         BlockedTensor& C1) {
    local_timer timer;

    if (C1.is_block("ac")) {
        auto C1ac = C1.block("ac");
        C1ac("um") += alpha * H1.block("aa")("vu") * T1.block("ca")("mv");
        C1ac("um") -= alpha * H1.block("cc")("mn") * T1.block("ca")("nu");
    }

    if (C1.is_block("vc")) {
        auto C1vc = C1.block("vc");
        C1vc("em") += alpha * H1.block("vv")("fe") * T1.block("cv")("mf");
        C1vc("em") -= alpha * H1.block("cc")("mn") * T1.block("cv")("ne");
    }

    if (C1.is_block("va")) {
        auto C1va = C1.block("va");
        C1va("eu") += alpha * H1.block("vv")("fe") * T1.block("av")("uf");
        C1va("eu") -= alpha * H1.block("aa")("uv") * T1.block("av")("ve");
    }

    if (C1.is_block("aa") and T1.is_block("aa")) {
        auto C1aa = C1.block("aa");
        auto H1aa = H1.block("aa");
        auto T1aa = T1.block("aa");
        C1aa("uv") += alpha * H1aa("xv") * T1aa("ux");
        C1aa("uv") += alpha * H1aa("ux") * T1aa("xv");
        C1aa("vu") += alpha * H1aa("xv") * T1aa("ux");
        C1aa("vu") += alpha * H1aa("ux") * T1aa("xv");
    }

    dsrg_time_.add("121", timer.get());
}

void SADSRG::H1d_A2_C1ph(BlockedTensor& H1, BlockedTensor& S2, const double& alpha,
                         BlockedTensor& C1) {
    /**
     * C1ph = [H1^{d}, A2]^{od}
     * C1["ia"] += 0.5 * alpha * H1["xu"] * S2["ivax"] * L1_["uv"];
     * C1["ia"] -= 0.5 * alpha * H1["vy"] * S2["iyau"] * L1_["uv"];
     *
     * S2["ijab"] = 2 * T2["ijab"] - T2["jiab"];
     * S2 blocks are assumed available in memory: aavv, ccaa, caav, acav, aava, caaa, aaaa
     */
    local_timer timer;
    auto H1d = H1.block("aa");
    auto na = actv_mos_.size();
    auto temp = ambit::Tensor::build(tensor_type_, "temp_H1d_A2_C1ph", {na, na});

    bool is_ac = C1.is_block("ac");
    bool is_vc = C1.is_block("vc");
    bool is_va = C1.is_block("va");
    bool is_aa = C1.is_block("aa");

    temp("xv") = H1d("xu") * L1_.block("aa")("uv");
    if (is_ac)
        C1.block("ac")("ym") += 0.5 * alpha * temp("xv") * S2.block("caaa")("mvyx");
    if (is_vc)
        C1.block("vc")("em") += 0.5 * alpha * temp("xv") * S2.block("acav")("vmxe");
    if (is_va)
        C1.block("va")("eu") += 0.5 * alpha * temp("xv") * S2.block("aava")("uvex");
    if (is_aa) {
        C1.block("aa")("uy") += 0.5 * alpha * temp("xv") * S2.block("aaaa")("uvyx");
        C1.block("aa")("yu") += 0.5 * alpha * temp("xv") * S2.block("aaaa")("uvyx");
    }

    temp("uy") = H1d("vy") * L1_.block("aa")("uv");
    if (is_ac)
        C1.block("ac")("vm") -= 0.5 * alpha * temp("uy") * S2.block("caaa")("myvu");
    if (is_va)
        C1.block("va")("ex") -= 0.5 * alpha * temp("uy") * S2.block("aava")("xyeu");
    if (is_vc)
        C1.block("vc")("em") -= 0.5 * alpha * temp("uy") * S2.block("acav")("ymue");
    if (is_aa) {
        C1.block("aa")("xv") -= 0.5 * alpha * temp("uy") * S2.block("aaaa")("xyvu");
        C1.block("aa")("vx") -= 0.5 * alpha * temp("uy") * S2.block("aaaa")("xyvu");
    }

    dsrg_time_.add("121", timer.get());
}

void SADSRG::H1d_A2_C2pphh(BlockedTensor& H1, BlockedTensor& T2, const double& alpha,
                           BlockedTensor& C2) {
    /**
     * The following blocks are assumed available in memory:
     * H1: cc, aa, vv
     * T2: hhpp including aaaa
     * C2: pphh including aaaa
     */
    auto H1c = H1.block("cc");
    auto H1a = H1.block("aa");
    auto H1v = H1.block("vv");
    auto nc = core_mos_.size();
    auto na = actv_mos_.size();
    auto nv = virt_mos_.size();

    // small T2 blocks
    std::vector<std::string> small_blocks;
    for (const auto block : {"vvaa", "aacc", "avca", "avac", "vaaa", "aaca", "aaaa"}) {
        std::string sblock;
        for (int i : {1, 0, 3, 2})
            sblock.push_back(block[i]);
        if (C2.is_block(block) or C2.is_block(sblock))
            small_blocks.push_back(block);
    }
    auto temp2 = ambit::BlockedTensor::build(tensor_type_, "temp_H1dA2C2pphh", small_blocks);
    H1d_A2_C2pphh_small(H1, T2, 1.0, temp2);
    C2["pqrs"] += alpha * temp2["pqrs"];
    for (const auto block : {"avca", "avac", "vaaa", "aaca"}) {
        std::string sblock;
        for (int i : {1, 0, 3, 2})
            sblock.push_back(block[i]);
        if (C2.is_block(sblock))
            C2.block(sblock)("qpsr") += alpha * temp2.block(block)("pqrs");
    }

    local_timer timer;

    // T2 ccvv
    if (C2.is_block("vvcc") and T2.is_block("ccvv")) {
        auto T2ccvv = T2.block("ccvv");
        auto temp2 = ambit::Tensor::build(tensor_type_, "temp_H1dA2C2_ccvv", {nc, nc, nv, nv});
        temp2("mnef") += H1v("ae") * T2ccvv("mnaf");
        temp2("mnef") -= H1c("mi") * T2ccvv("inef");
        C2.block("vvcc")("efmn") += alpha * temp2("mnef");
        C2.block("vvcc")("fenm") += alpha * temp2("mnef");
    }

    // T2 ccva or ccav
    if (C2.is_block("vacc") or C2.is_block("avcc")) {
        auto temp2 = ambit::Tensor::build(tensor_type_, "temp_H1dA2C2_ccav", {nc, nc, na, nv});
        if (T2.is_block("ccav")) {
            auto T2ccav = T2.block("ccav");
            temp2("mnue") += H1a("au") * T2ccav("mnae");
            temp2("mnue") += H1v("ae") * T2ccav("mnua");
            temp2("mnue") -= H1c("mi") * T2ccav("inue");
            temp2("mnue") -= H1c("ni") * T2ccav("miue");
        } else if (T2.is_block("ccva")) {
            auto T2ccva = T2.block("ccva");
            temp2("mnue") += H1a("au") * T2ccva("nmea");
            temp2("mnue") += H1v("ae") * T2ccva("nmau");
            temp2("mnue") -= H1c("mi") * T2ccva("nieu");
            temp2("mnue") -= H1c("ni") * T2ccva("imeu");
        }
        if (C2.is_block("vacc"))
            C2.block("vacc")("eunm") += alpha * temp2("mnue");
        if (C2.is_block("avcc"))
            C2.block("avcc")("uemn") += alpha * temp2("mnue");
    }

    // T2 acvv or cavv
    if (C2.is_block("vvca") or C2.is_block("vvac")) {
        auto temp2 = ambit::Tensor::build(tensor_type_, "temp_H1dA2C2_acvv", {na, nc, nv, nv});
        if (T2.is_block("acvv")) {
            auto T2acvv = T2.block("acvv");
            temp2("umef") -= H1a("ui") * T2acvv("imef");
            temp2("umef") -= H1c("mi") * T2acvv("uief");
            temp2("umef") += H1v("ae") * T2acvv("umaf");
            temp2("umef") += H1v("af") * T2acvv("umea");
        } else if (T2.is_block("cavv")) {
            auto T2cavv = T2.block("cavv");
            temp2("umef") -= H1a("ui") * T2cavv("mife");
            temp2("umef") -= H1c("mi") * T2cavv("iufe");
            temp2("umef") += H1v("ae") * T2cavv("mufa");
            temp2("umef") += H1v("af") * T2cavv("muae");
        }
        if (C2.is_block("vvac"))
            C2.block("vvac")("efum") += alpha * temp2("umef");
        if (C2.is_block("vvca"))
            C2.block("vvca")("femu") += alpha * temp2("umef");
    }

    dsrg_time_.add("122", timer.get());
}

void SADSRG::H1d_A2_C2pphh_small(BlockedTensor& H1, BlockedTensor& T2, const double& alpha,
                                 BlockedTensor& C2) {
    /**
     * C2["ijpb"] += alpha * T2["ijab"] * H1["ap"];
     * C2["jibp"] += alpha * T2["ijab"] * H1["ap"];
     * C2["qjab"] -= alpha * T2["ijab"] * H1["qi"];
     * C2["jqba"] -= alpha * T2["ijab"] * H1["qi"];
     *
     * The following blocks should be available in memory:
     * H1: cc, aa, vv
     * T2: aavv, ccaa, caav, acav, aava, caaa, aaaa
     * C2: vvaa, aacc, avca, avac, vaaa, aaca, aaaa
     */
    local_timer timer;

    auto H1c = H1.block("cc");
    auto H1a = H1.block("aa");
    auto H1v = H1.block("vv");

    auto nc = core_mos_.size();
    auto na = actv_mos_.size();
    auto nv = virt_mos_.size();

    ambit::Tensor temp;
    std::string prefix = "temp_H1dA2C2pphh";

    if (C2.is_block("vvaa") and T2.is_block("aavv")) {
        temp = ambit::Tensor::build(tensor_type_, prefix + "_aavv", {na, na, nv, nv});
        temp("uvef") += T2.block("aavv")("uvaf") * H1v("ae");
        temp("uvef") -= T2.block("aavv")("ivef") * H1a("ui");
        C2.block("vvaa")("efuv") += alpha * temp("uvef");
        C2.block("vvaa")("fevu") += alpha * temp("uvef");
    }

    if (C2.is_block("aacc") and T2.is_block("ccaa")) {
        temp = ambit::Tensor::build(tensor_type_, prefix + "_ccaa", {nc, nc, na, na});
        temp("mnuv") -= T2.block("ccaa")("inuv") * H1c("mi");
        temp("mnuv") += T2.block("ccaa")("mnav") * H1a("au");
        C2.block("aacc")("uvmn") += alpha * temp("mnuv");
        C2.block("aacc")("vunm") += alpha * temp("mnuv");
    }

    if (C2.is_block("avca") and T2.is_block("caav")) {
        auto C2avca = C2.block("avca");
        auto T2caav = T2.block("caav");
        C2avca("vemu") += alpha * T2caav("muvf") * H1v("fe");
        C2avca("vemu") += alpha * T2caav("muye") * H1a("yv");
        C2avca("vemu") -= alpha * T2caav("nuve") * H1c("mn");
        C2avca("vemu") -= alpha * T2caav("mxve") * H1a("ux");
    }

    if (C2.is_block("avac") and T2.is_block("acav")) {
        auto C2avac = C2.block("avac");
        auto T2acav = T2.block("acav");
        C2avac("veum") += alpha * T2acav("umvf") * H1v("fe");
        C2avac("veum") += alpha * T2acav("umye") * H1a("yv");
        C2avac("veum") -= alpha * T2acav("unve") * H1c("mn");
        C2avac("veum") -= alpha * T2acav("xmve") * H1a("ux");
    }

    if (C2.is_block("vaaa") and T2.is_block("aava")) {
        auto C2vaaa = C2.block("vaaa");
        auto T2aava = T2.block("aava");
        C2vaaa("eyuv") += alpha * T2aava("uvfy") * H1v("fe");
        C2vaaa("eyuv") += alpha * T2aava("uvex") * H1a("xy");
        C2vaaa("eyuv") -= alpha * T2aava("xvey") * H1a("ux");
        C2vaaa("eyuv") -= alpha * T2aava("uxey") * H1a("vx");
    }

    if (C2.is_block("aaca") and T2.is_block("caaa")) {
        auto C2aaca = C2.block("aaca");
        auto T2caaa = T2.block("caaa");
        C2aaca("xymu") -= alpha * T2caaa("nuxy") * H1c("mn");
        C2aaca("xymu") -= alpha * T2caaa("mvxy") * H1a("uv");
        C2aaca("xymu") += alpha * T2caaa("muvy") * H1a("vx");
        C2aaca("xymu") += alpha * T2caaa("muxv") * H1a("vy");
    }

    if (C2.is_block("aaaa") and T2.is_block("aaaa")) {
        temp = ambit::Tensor::build(tensor_type_, prefix + "_aaaa", {na, na, na, na});
        temp("uvxy") += T2.block("aaaa")("uvay") * H1a("ax");
        temp("uvxy") -= T2.block("aaaa")("ivxy") * H1a("ui");
        C2.block("aaaa")("xyuv") += alpha * temp("uvxy");
        C2.block("aaaa")("yxvu") += alpha * temp("uvxy");
        C2.block("aaaa")("uvxy") += alpha * temp("uvxy");
        C2.block("aaaa")("vuyx") += alpha * temp("uvxy");
    }

    dsrg_time_.add("122", timer.get());
}

void SADSRG::H1_A1_C1ph(BlockedTensor& H1, BlockedTensor& T1, const double& alpha,
                        BlockedTensor& C1) {
    local_timer timer;

    C1["ai"] += alpha * H1["ca"] * T1["ic"];
    C1["ai"] -= alpha * H1["ik"] * T1["ka"];

    C1["ui"] += alpha * H1["ai"] * T1["ua"];
    C1["au"] -= alpha * H1["ai"] * T1["iu"];

    dsrg_time_.add("111", timer.get());
}

void SADSRG::H1_A2_C1ph(BlockedTensor& H1, BlockedTensor& S2, const double& alpha,
                        BlockedTensor& C1) {
    local_timer timer;

    C1["ai"] += alpha * H1["bm"] * S2["imab"];
    C1["ai"] += 0.5 * alpha * H1["bu"] * S2["ivab"] * L1_["uv"];
    C1["ai"] -= 0.5 * alpha * H1["vj"] * S2["ijau"] * L1_["uv"];

    if (C1.is_block("aa")) {
        C1["xy"] += alpha * H1["bm"] * S2["xmyb"];
        C1["xy"] += 0.5 * alpha * H1["bu"] * S2["xvyb"] * L1_["uv"];
        C1["xy"] -= 0.5 * alpha * H1["vj"] * S2["xjyu"] * L1_["uv"];
    }

    dsrg_time_.add("121", timer.get());
}
} // namespace forte
