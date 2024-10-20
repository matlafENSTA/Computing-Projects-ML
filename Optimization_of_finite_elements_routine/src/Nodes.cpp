#include "Nodes.hpp"

// ======================================================================
// Build local coordinates of nodes on the reference element
//   Input: N
//   Output: _r, _s, _t
// ======================================================================

void getReferenceNodes(int N, vector<double> &_r, vector<double> &_s, vector<double> &_t) {

  cout << "CALL getReferenceNodes()\n";

  double ref_r_01[4] = {-1, 1, -1, -1};
  double ref_s_01[4] = {-1, -1, 1, -1};
  double ref_t_01[4] = {-1, -1, -1, 1};

  double ref_r_02[10] = {-1, 0, 1, -1, 0, -1, -1, 0, -1, -1};
  double ref_s_02[10] = {-1, -1, -1, 0, 0, 1, -1, -1, 0, -1};
  double ref_t_02[10] = {-1, -1, -1, -1, -1, -1, 0, 0, 0, 1};

  double ref_r_03[20] = {-1, -0.447213595499958, 0.447213595499958, 1, -1, -0.333333333333333, 0.447213595499958, -1, -0.447213595499958, -1, -1, -0.333333333333333, 0.447213595499958, -1, -0.333333333333333, -1, -1, -0.447213595499958, -1, -1};
  double ref_s_03[20] = {-1, -1, -1, -1, -0.447213595499958, -0.333333333333333, -0.447213595499958, 0.447213595499958, 0.447213595499958, 1, -1, -1, -1, -0.333333333333333, -0.333333333333333, 0.447213595499958, -1, -1, -0.447213595499958, -1};
  double ref_t_03[20] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -0.447213595499958, -0.333333333333333, -0.447213595499958, -0.333333333333333, -0.333333333333333, -0.447213595499958, 0.447213595499958, 0.447213595499958, 0.447213595499958, 1};

  double ref_r_04[35] = {-1, -0.654653670707977, 0, 0.654653670707977, 1, -1, -0.551583572090994, 0.103167144181987, 0.654653670707977, -1, -0.551583572090994, 0, -1, -0.654653670707977, -1, -1, -0.551583572090994, 0.103167144181987, 0.654653670707977, -1, -0.5, 0.103167144181987, -1, -0.551583572090994, -1, -1, -0.551583572090994, 0, -1, -0.551583572090994, -1, -1, -0.654653670707977, -1, -1};
  double ref_s_04[35] = {-1, -1, -1, -1, -1, -0.654653670707977, -0.551583572090994, -0.551583572090994, -0.654653670707977, 0, 0.103167144181987, 0, 0.654653670707977, 0.654653670707977, 1, -1, -1, -1, -1, -0.551583572090994, -0.5, -0.551583572090994, 0.103167144181987, 0.103167144181987, 0.654653670707977, -1, -1, -1, -0.551583572090993, -0.551583572090993, 0, -1, -1, -0.654653670707977, -1};
  double ref_t_04[35] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -0.654653670707977, -0.551583572090993, -0.551583572090994, -0.654653670707977, -0.551583572090993, -0.5, -0.551583572090993, -0.551583572090994, -0.551583572090994, -0.654653670707977, 0, 0.103167144181987, 0, 0.103167144181987, 0.103167144181987, 0, 0.654653670707977, 0.654653670707977, 0.654653670707977, 1};

  double ref_r_05[56] = {-1, -0.765055323929465, -0.285231516480645, 0.285231516480645, 0.765055323929465, 1, -1, -0.68854346464994, -0.165752193680209, 0.377086929299879, 0.765055323929465, -1, -0.668495612639582, -0.165752193680209, 0.285231516480645, -1, -0.68854346464994, -0.285231516480645, -1, -0.765055323929465, -1, -1, -0.68854346464994, -0.165752193680209, 0.377086929299879, 0.765055323929465, -1, -0.622331815219382, -0.133004554341855, 0.377086929299879, -1, -0.622331815219382, -0.165752193680209, -1, -0.68854346464994, -1, -1, -0.668495612639582, -0.165752193680209, 0.285231516480645, -1, -0.622331815219382, -0.165752193680209, -1, -0.668495612639582, -1, -1, -0.68854346464994, -0.285231516480645, -1, -0.68854346464994, -1, -1, -0.765055323929465, -1, -1};
  double ref_s_05[56] = {-1, -1, -1, -1, -1, -1, -0.765055323929465, -0.68854346464994, -0.668495612639582, -0.68854346464994, -0.765055323929465, -0.285231516480645, -0.165752193680209, -0.165752193680209, -0.285231516480645, 0.285231516480645, 0.377086929299879, 0.285231516480645, 0.765055323929465, 0.765055323929465, 1, -1, -1, -1, -1, -1, -0.68854346464994, -0.622331815219382, -0.622331815219382, -0.68854346464994, -0.165752193680209, -0.133004554341855, -0.165752193680209, 0.377086929299879, 0.377086929299879, 0.765055323929465, -1, -1, -1, -1, -0.668495612639582, -0.622331815219382, -0.668495612639582, -0.165752193680209, -0.165752193680209, 0.285231516480645, -1, -1, -1, -0.688543464649939, -0.688543464649939, -0.285231516480645, -1, -1, -0.765055323929465, -1};
  double ref_t_05[56] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -0.765055323929465, -0.688543464649939, -0.668495612639582, -0.688543464649939, -0.765055323929465, -0.688543464649939, -0.622331815219382, -0.622331815219382, -0.688543464649939, -0.668495612639582, -0.622331815219382, -0.668495612639582, -0.688543464649939, -0.688543464649939, -0.765055323929465, -0.285231516480645, -0.165752193680209, -0.165752193680209, -0.285231516480645, -0.165752193680209, -0.133004554341855, -0.165752193680209, -0.165752193680209, -0.165752193680209, -0.285231516480645, 0.285231516480645, 0.377086929299879, 0.285231516480645, 0.377086929299879, 0.377086929299879, 0.285231516480645, 0.765055323929465, 0.765055323929465, 0.765055323929465, 1};

  double ref_r_06[84] = {-1, -0.830223896278567, -0.468848793470714, 0, 0.468848793470714, 0.830223896278567, 1, -1, -0.773622922202285, -0.360567121835569, 0.119298207804603, 0.547245844404568, 0.830223896278567, -1, -0.758731085969035, -0.333333333333333, 0.119298207804604, 0.468848793470714, -1, -0.758731085969035, -0.360567121835569, 8.58760498152762e-18, -1, -0.773622922202284, -0.468848793470714, -1, -0.830223896278567, -1, -1, -0.773622922202285, -0.360567121835569, 0.119298207804603, 0.547245844404568, 0.830223896278567, -1, -0.710802806967777, -0.305827490438019, 0.132408420903329, 0.547245844404568, -1, -0.694172509561981, -0.305827490438019, 0.119298207804603, -1, -0.710802806967777, -0.360567121835569, -1, -0.773622922202284, -1, -1, -0.758731085969035, -0.333333333333333, 0.119298207804604, 0.468848793470714, -1, -0.694172509561981, -0.305827490438019, 0.119298207804603, -1, -0.694172509561981, -0.333333333333333, -1, -0.758731085969035, -1, -1, -0.758731085969035, -0.360567121835569, 0, -1, -0.710802806967777, -0.360567121835569, -1, -0.758731085969035, -1, -1, -0.773622922202284, -0.468848793470714, -1, -0.773622922202284, -1, -1, -0.830223896278567, -1, -1};
  double ref_s_06[84] = {-1, -1, -1, -1, -1, -1, -1, -0.830223896278567, -0.773622922202284, -0.758731085969035, -0.758731085969035, -0.773622922202284, -0.830223896278567, -0.468848793470714, -0.360567121835569, -0.333333333333333, -0.360567121835569, -0.468848793470714, -1.28197512425571e-16, 0.119298207804603, 0.119298207804603, -1.28197512425571e-16, 0.468848793470714, 0.547245844404568, 0.468848793470714, 0.830223896278567, 0.830223896278567, 1, -1, -1, -1, -1, -1, -1, -0.773622922202284, -0.710802806967777, -0.694172509561981, -0.710802806967777, -0.773622922202284, -0.360567121835569, -0.305827490438019, -0.305827490438019, -0.360567121835569, 0.119298207804603, 0.132408420903329, 0.119298207804603, 0.547245844404568, 0.547245844404568, 0.830223896278567, -1, -1, -1, -1, -1, -0.758731085969035, -0.694172509561981, -0.694172509561981, -0.758731085969035, -0.333333333333333, -0.305827490438019, -0.333333333333333, 0.119298207804604, 0.119298207804604, 0.468848793470714, -1, -1, -1, -1, -0.758731085969035, -0.710802806967776, -0.758731085969035, -0.360567121835569, -0.360567121835569, 2.2662332591842e-17, -1, -1, -1, -0.773622922202284, -0.773622922202284, -0.468848793470714, -1, -1, -0.830223896278567, -1};
  double ref_t_06[84] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -0.830223896278567, -0.773622922202284, -0.758731085969035, -0.758731085969035, -0.773622922202284, -0.830223896278567, -0.773622922202284, -0.710802806967776, -0.694172509561981, -0.710802806967776, -0.773622922202284, -0.758731085969035, -0.694172509561981, -0.694172509561981, -0.758731085969035, -0.758731085969035, -0.710802806967776, -0.758731085969035, -0.773622922202284, -0.773622922202284, -0.830223896278567, -0.468848793470714, -0.360567121835569, -0.333333333333333, -0.360567121835569, -0.468848793470714, -0.360567121835569, -0.305827490438019, -0.305827490438019, -0.360567121835569, -0.333333333333333, -0.305827490438019, -0.333333333333333, -0.360567121835569, -0.360567121835569, -0.468848793470714, 6.79869977755259e-17, 0.119298207804604, 0.119298207804604, 6.79869977755259e-17, 0.119298207804604, 0.132408420903329, 0.119298207804604, 0.119298207804604, 0.119298207804603, -6.79869977755259e-17, 0.468848793470714, 0.547245844404568, 0.468848793470714, 0.547245844404568, 0.547245844404568, 0.468848793470714, 0.830223896278567, 0.830223896278567, 0.830223896278567, 1};

  double ref_r_07[120] = {-1, -0.871740148509607, -0.591700181433142, -0.209299217902479, 0.209299217902479, 0.591700181433142, 0.871740148509607, 1, -1, -0.825283816702066, -0.502211370361116, -0.0989582963347835, 0.30891067581573, 0.650567633404132, 0.871740148509606, -1, -0.806699305454614, -0.466671445131973, -0.0666571097360542, 0.30891067581573, 0.591700181433142, -1, -0.802083407330434, -0.466671445131973, -0.0989582963347831, 0.209299217902479, -1, -0.806699305454614, -0.502211370361116, -0.209299217902479, -1, -0.825283816702066, -0.591700181433142, -1, -0.871740148509607, -1, -1, -0.825283816702066, -0.502211370361116, -0.0989582963347835, 0.30891067581573, 0.650567633404132, 0.871740148509607, -1, -0.769697876091516, -0.444458691747615, -0.0596324566484163, 0.309093628274547, 0.650567633404132, -1, -0.747954425801985, -0.42209801484566, -0.059632456648416, 0.30891067581573, -1, -0.747954425801985, -0.444458691747615, -0.0989582963347832, -1, -0.769697876091516, -0.502211370361116, -1, -0.825283816702066, -1, -1, -0.806699305454614, -0.466671445131973, -0.0666571097360542, 0.30891067581573, 0.591700181433142, -1, -0.747954425801985, -0.42209801484566, -0.0596324566484162, 0.30891067581573, -1, -0.733705955463021, -0.42209801484566, -0.0666571097360542, -1, -0.747954425801985, -0.466671445131973, -1, -0.806699305454615, -1, -1, -0.802083407330434, -0.466671445131973, -0.0989582963347831, 0.209299217902479, -1, -0.747954425801985, -0.444458691747615, -0.0989582963347834, -1, -0.747954425801985, -0.466671445131973, -1, -0.802083407330434, -1, -1, -0.806699305454615, -0.502211370361116, -0.209299217902479, -1, -0.769697876091516, -0.502211370361116, -1, -0.806699305454615, -1, -1, -0.825283816702066, -0.591700181433142, -1, -0.825283816702066, -1, -1, -0.871740148509607, -1, -1};
  double ref_s_07[120] = {-1, -1, -1, -1, -1, -1, -1, -1, -0.871740148509607, -0.825283816702066, -0.806699305454614, -0.802083407330433, -0.806699305454614, -0.825283816702066, -0.871740148509606, -0.591700181433142, -0.502211370361116, -0.466671445131973, -0.466671445131973, -0.502211370361116, -0.591700181433142, -0.209299217902479, -0.0989582963347832, -0.0666571097360543, -0.0989582963347834, -0.209299217902479, 0.209299217902479, 0.30891067581573, 0.30891067581573, 0.209299217902479, 0.591700181433142, 0.650567633404131, 0.591700181433142, 0.871740148509606, 0.871740148509606, 1, -1, -1, -1, -1, -1, -1, -1, -0.825283816702066, -0.769697876091516, -0.747954425801985, -0.747954425801985, -0.769697876091516, -0.825283816702066, -0.502211370361116, -0.444458691747615, -0.42209801484566, -0.444458691747615, -0.502211370361116, -0.0989582963347835, -0.0596324566484162, -0.0596324566484162, -0.0989582963347835, 0.30891067581573, 0.309093628274547, 0.30891067581573, 0.650567633404132, 0.650567633404132, 0.871740148509606, -1, -1, -1, -1, -1, -1, -0.806699305454614, -0.747954425801985, -0.733705955463021, -0.747954425801985, -0.806699305454614, -0.466671445131973, -0.42209801484566, -0.42209801484566, -0.466671445131973, -0.0666571097360542, -0.0596324566484161, -0.0666571097360542, 0.30891067581573, 0.30891067581573, 0.591700181433142, -1, -1, -1, -1, -1, -0.802083407330433, -0.747954425801985, -0.747954425801985, -0.802083407330433, -0.466671445131973, -0.444458691747614, -0.466671445131973, -0.0989582963347831, -0.0989582963347831, 0.209299217902479, -1, -1, -1, -1, -0.806699305454614, -0.769697876091516, -0.806699305454614, -0.502211370361116, -0.502211370361116, -0.209299217902479, -1, -1, -1, -0.825283816702066, -0.825283816702066, -0.591700181433142, -1, -1, -0.871740148509606, -1};
  double ref_t_07[120] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -0.871740148509607, -0.825283816702066, -0.806699305454614, -0.802083407330433, -0.806699305454614, -0.825283816702066, -0.871740148509606, -0.825283816702066, -0.769697876091516, -0.747954425801985, -0.747954425801985, -0.769697876091516, -0.825283816702066, -0.806699305454614, -0.747954425801985, -0.733705955463021, -0.747954425801985, -0.806699305454615, -0.802083407330433, -0.747954425801985, -0.747954425801985, -0.802083407330433, -0.806699305454614, -0.769697876091516, -0.806699305454614, -0.825283816702066, -0.825283816702066, -0.871740148509606, -0.591700181433142, -0.502211370361116, -0.466671445131973, -0.466671445131973, -0.502211370361115, -0.591700181433142, -0.502211370361116, -0.444458691747614, -0.42209801484566, -0.444458691747614, -0.502211370361116, -0.466671445131973, -0.42209801484566, -0.42209801484566, -0.466671445131973, -0.466671445131973, -0.444458691747614, -0.466671445131973, -0.502211370361115, -0.502211370361115, -0.591700181433142, -0.209299217902479, -0.0989582963347832, -0.0666571097360542, -0.0989582963347834, -0.209299217902479, -0.0989582963347832, -0.0596324566484161, -0.0596324566484162, -0.0989582963347832, -0.0666571097360542, -0.0596324566484162, -0.0666571097360542, -0.0989582963347834, -0.0989582963347834, -0.209299217902479, 0.209299217902479, 0.30891067581573, 0.30891067581573, 0.209299217902479, 0.30891067581573, 0.309093628274547, 0.30891067581573, 0.30891067581573, 0.30891067581573, 0.209299217902479, 0.591700181433142, 0.650567633404132, 0.591700181433142, 0.650567633404132, 0.650567633404132, 0.591700181433143, 0.871740148509607, 0.871740148509606, 0.871740148509607, 1};

  double ref_r_08[165] = {-1, -0.89975799541146, -0.677186279510738, -0.363117463826178, 0, 0.363117463826178, 0.677186279510738, 0.89975799541146, 1, -1, -0.861709494494173, -0.602621886086653, -0.265044207129911, 0.10143688572376, 0.446014634733998, 0.723418988988346, 0.89975799541146, -1, -0.843392748647346, -0.566293655878678, -0.222195622530102, 0.132587311757356, 0.446014634733998, 0.677186279510738, -1, -0.836392678593849, -0.555608754939795, -0.222195622530103, 0.10143688572376, 0.363117463826178, -1, -0.836392678593849, -0.566293655878678, -0.265044207129911, 0, -1, -0.843392748647345, -0.602621886086653, -0.363117463826178, -1, -0.861709494494173, -0.677186279510738, -1, -0.89975799541146, -1, -1, -0.861709494494173, -0.602621886086653, -0.265044207129911, 0.10143688572376, 0.446014634733998, 0.723418988988346, 0.89975799541146, -1, -0.812644144069609, -0.54857984257542, -0.216193164886725, 0.128143146183281, 0.437932432208828, 0.723418988988346, -1, -0.78978165180393, -0.51657615948891, -0.195436130795035, 0.128143146183281, 0.446014634733998, -1, -0.783806835113275, -0.516576159488909, -0.216193164886725, 0.10143688572376, -1, -0.78978165180393, -0.54857984257542, -0.265044207129911, -1, -0.812644144069609, -0.602621886086653, -1, -0.861709494494173, -1, -1, -0.843392748647346, -0.566293655878678, -0.222195622530102, 0.132587311757356, 0.446014634733998, 0.677186279510738, -1, -0.78978165180393, -0.51657615948891, -0.195436130795035, 0.128143146183281, 0.446014634733998, -1, -0.771411550227146, -0.5, -0.195436130795035, 0.132587311757356, -1, -0.771411550227146, -0.516576159488909, -0.222195622530102, -1, -0.78978165180393, -0.566293655878678, -1, -0.843392748647345, -1, -1, -0.836392678593849, -0.555608754939795, -0.222195622530103, 0.10143688572376, 0.363117463826178, -1, -0.783806835113275, -0.516576159488909, -0.216193164886725, 0.10143688572376, -1, -0.771411550227146, -0.51657615948891, -0.222195622530103, -1, -0.783806835113275, -0.555608754939795, -1, -0.836392678593849, -1, -1, -0.836392678593849, -0.566293655878678, -0.265044207129911, 0, -1, -0.78978165180393, -0.54857984257542, -0.265044207129911, -1, -0.789781651803931, -0.566293655878678, -1, -0.836392678593849, -1, -1, -0.843392748647345, -0.602621886086653, -0.363117463826178, -1, -0.812644144069609, -0.602621886086653, -1, -0.843392748647346, -1, -1, -0.861709494494173, -0.677186279510738, -1, -0.861709494494173, -1, -1, -0.89975799541146, -1, -1};
  double ref_s_08[165] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -0.89975799541146, -0.861709494494173, -0.843392748647345, -0.836392678593849, -0.836392678593849, -0.843392748647346, -0.861709494494173, -0.89975799541146, -0.677186279510738, -0.602621886086653, -0.566293655878678, -0.555608754939795, -0.566293655878678, -0.602621886086653, -0.677186279510738, -0.363117463826178, -0.265044207129911, -0.222195622530102, -0.222195622530102, -0.265044207129911, -0.363117463826178, 0, 0.10143688572376, 0.132587311757356, 0.10143688572376, 0, 0.363117463826178, 0.446014634733998, 0.446014634733998, 0.363117463826178, 0.677186279510738, 0.723418988988346, 0.677186279510738, 0.89975799541146, 0.89975799541146, 1, -1, -1, -1, -1, -1, -1, -1, -1, -0.861709494494173, -0.812644144069609, -0.789781651803931, -0.783806835113275, -0.789781651803931, -0.812644144069609, -0.861709494494173, -0.602621886086653, -0.54857984257542, -0.516576159488909, -0.516576159488909, -0.54857984257542, -0.602621886086653, -0.265044207129911, -0.216193164886725, -0.195436130795035, -0.216193164886725, -0.265044207129911, 0.10143688572376, 0.128143146183281, 0.128143146183281, 0.10143688572376, 0.446014634733998, 0.437932432208828, 0.446014634733998, 0.723418988988346, 0.723418988988346, 0.89975799541146, -1, -1, -1, -1, -1, -1, -1, -0.843392748647345, -0.789781651803931, -0.771411550227146, -0.771411550227146, -0.78978165180393, -0.843392748647345, -0.566293655878678, -0.516576159488909, -0.5, -0.516576159488909, -0.566293655878678, -0.222195622530102, -0.195436130795035, -0.195436130795035, -0.222195622530102, 0.132587311757356, 0.128143146183281, 0.132587311757356, 0.446014634733998, 0.446014634733998, 0.677186279510738, -1, -1, -1, -1, -1, -1, -0.836392678593849, -0.783806835113275, -0.771411550227146, -0.783806835113275, -0.836392678593849, -0.555608754939795, -0.516576159488909, -0.516576159488909, -0.555608754939795, -0.222195622530102, -0.216193164886725, -0.222195622530102, 0.10143688572376, 0.10143688572376, 0.363117463826178, -1, -1, -1, -1, -1, -0.836392678593849, -0.78978165180393, -0.78978165180393, -0.836392678593849, -0.566293655878678, -0.54857984257542, -0.566293655878678, -0.265044207129911, -0.265044207129911, 0, -1, -1, -1, -1, -0.843392748647345, -0.812644144069609, -0.843392748647345, -0.602621886086653, -0.602621886086653, -0.363117463826178, -1, -1, -1, -0.861709494494173, -0.861709494494173, -0.677186279510738, -1, -1, -0.89975799541146, -1};
  double ref_t_08[165] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -0.89975799541146, -0.861709494494173, -0.843392748647345, -0.836392678593849, -0.836392678593849, -0.843392748647345, -0.861709494494173, -0.89975799541146, -0.861709494494173, -0.812644144069609, -0.78978165180393, -0.783806835113275, -0.78978165180393, -0.812644144069609, -0.861709494494173, -0.843392748647345, -0.78978165180393, -0.771411550227146, -0.771411550227146, -0.78978165180393, -0.843392748647345, -0.836392678593849, -0.783806835113275, -0.771411550227146, -0.783806835113275, -0.836392678593849, -0.836392678593849, -0.78978165180393, -0.78978165180393, -0.836392678593849, -0.843392748647345, -0.812644144069609, -0.843392748647345, -0.861709494494173, -0.861709494494173, -0.89975799541146, -0.677186279510738, -0.602621886086653, -0.566293655878678, -0.555608754939795, -0.566293655878678, -0.602621886086653, -0.677186279510738, -0.602621886086653, -0.54857984257542, -0.516576159488909, -0.516576159488909, -0.54857984257542, -0.602621886086653, -0.566293655878678, -0.516576159488909, -0.5, -0.516576159488909, -0.566293655878678, -0.555608754939795, -0.516576159488909, -0.516576159488909, -0.555608754939795, -0.566293655878678, -0.54857984257542, -0.566293655878678, -0.602621886086653, -0.602621886086653, -0.677186279510738, -0.363117463826178, -0.265044207129911, -0.222195622530102, -0.222195622530102, -0.265044207129911, -0.363117463826178, -0.265044207129911, -0.216193164886725, -0.195436130795035, -0.216193164886725, -0.265044207129911, -0.222195622530102, -0.195436130795035, -0.195436130795035, -0.222195622530102, -0.222195622530102, -0.216193164886725, -0.222195622530102, -0.265044207129911, -0.265044207129911, -0.363117463826178, 0, 0.10143688572376, 0.132587311757356, 0.10143688572376, 0, 0.10143688572376, 0.128143146183281, 0.128143146183281, 0.10143688572376, 0.132587311757356, 0.128143146183281, 0.132587311757356, 0.10143688572376, 0.10143688572376, 0, 0.363117463826178, 0.446014634733998, 0.446014634733998, 0.363117463826178, 0.446014634733998, 0.437932432208828, 0.446014634733998, 0.446014634733998, 0.446014634733998, 0.363117463826178, 0.677186279510738, 0.723418988988346, 0.677186279510738, 0.723418988988346, 0.723418988988346, 0.677186279510738, 0.89975799541146, 0.89975799541146, 0.89975799541146, 1};

  double *ref_r;
  double *ref_s;
  double *ref_t;
  switch (N) {
  case 1:
    ref_r = ref_r_01;
    ref_s = ref_s_01;
    ref_t = ref_t_01;
    break;
  case 2:
    ref_r = ref_r_02;
    ref_s = ref_s_02;
    ref_t = ref_t_02;
    break;
  case 3:
    ref_r = ref_r_03;
    ref_s = ref_s_03;
    ref_t = ref_t_03;
    break;
  case 4:
    ref_r = ref_r_04;
    ref_s = ref_s_04;
    ref_t = ref_t_04;
    break;
  case 5:
    ref_r = ref_r_05;
    ref_s = ref_s_05;
    ref_t = ref_t_05;
    break;
  case 6:
    ref_r = ref_r_06;
    ref_s = ref_s_06;
    ref_t = ref_t_06;
    break;
  case 7:
    ref_r = ref_r_07;
    ref_s = ref_s_07;
    ref_t = ref_t_07;
    break;
  case 8:
    ref_r = ref_r_08;
    ref_s = ref_s_08;
    ref_t = ref_t_08;
    break;
  }

  int Np = (N + 1) * (N + 2) * (N + 3) / 6;
  _r.resize(Np);
  _s.resize(Np);
  _t.resize(Np);
  for (int n = 0; n < Np; n++) {
    _r[n] = ref_r[n];
    _s[n] = ref_s[n];
    _t[n] = ref_t[n];
  }
}

// ======================================================================
// Build global coordinates of nodes on the physical elements
//   Input: _EToV, _VX, _VY, _VZ, _r, _s, _t
//   Output: _x, _y, _z
// ======================================================================

void getPhysicalNodes(vector<int> &_EToV, vector<double> &_VX, vector<double> &_VY, vector<double> &_VZ, vector<double> &_r, vector<double> &_s, vector<double> &_t, vector<double> &_x, vector<double> &_y, vector<double> &_z) {

  cout << "CALL getPhysicalNodes()\n";

  int K = _EToV.size()/NVertTet;
  int Np = _r.size();

  _x.resize(K * Np);
  _y.resize(K * Np);
  _z.resize(K * Np);

  for (int k = 0, cntv = 0, cntf = 0; k < K; ++k) {
    int vGlo1 = _EToV[k * NVertTet + 0];
    int vGlo2 = _EToV[k * NVertTet + 1];
    int vGlo3 = _EToV[k * NVertTet + 2];
    int vGlo4 = _EToV[k * NVertTet + 3];
    double x1 = _VX[vGlo1], y1 = _VY[vGlo1], z1 = _VZ[vGlo1];
    double x2 = _VX[vGlo2], y2 = _VY[vGlo2], z2 = _VZ[vGlo2];
    double x3 = _VX[vGlo3], y3 = _VY[vGlo3], z3 = _VZ[vGlo3];
    double x4 = _VX[vGlo4], y4 = _VY[vGlo4], z4 = _VZ[vGlo4];
    for (int n = 0; n < Np; ++n) {
      double lam1 = -(1 + _r[n] + _s[n] + _t[n]);
      double lam2 = 1 + _r[n];
      double lam3 = 1 + _s[n];
      double lam4 = 1 + _t[n];
      _x[k * Np + n] = 0.5 * (lam1 * x1 + lam2 * x2 + lam3 * x3 + lam4 * x4);
      _y[k * Np + n] = 0.5 * (lam1 * y1 + lam2 * y2 + lam3 * y3 + lam4 * y4);
      _z[k * Np + n] = 0.5 * (lam1 * z1 + lam2 * z2 + lam3 * z3 + lam4 * z4);
      ++cntv;
    }
  }
}