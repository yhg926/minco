#ifndef MODEL_ANI_H
#define MODEL_ANI_H

#include <math.h>
#include "sketch_rearrange.h"

// model 1: 3-way linear model parameters

// optimized parameters for 3 denominator of 3-way linear model
// NUM_CODENS==11 values are recalibrated for minco 10k context-minhash
// sketches against Nayfach ANIm labels.
static const double opt_denom_params4x3[12] =
    {
#if (NUM_CODENS == 9 || NUM_CODENS == 10)
        1.0090536851341, 1.0013604061893, 1.03540226374509, 1.04270754727088, // denom1
        1.0108940845453, 1.0010797461052, 0.98914550652945, 1.01536320288592, // denom2
        1.0061957741437, 0.9957814716168, 0.98352379677340, 1.00951820561544  // denom3
#elif (NUM_CODENS == 11)
        1, 0.131595734232788, 1.54251229924664, 0.116283436444129,    // denom1
        1, 0.456844692517648, 1.44359887460418, 0.493945860006579,    // denom2
        1, 1.89001465044101, -1.98206480156227, -1.98562980351413     // denom3
#endif
};
// ANI>95 subset parameters
static const double N95opt_denom_params4x3[12] =
    {
#if (NUM_CODENS == 9 || NUM_CODENS == 10)
        1.0196821174688, 0.9817771341602, 1.04637038371202, 1.02441485561367, // denom1
        1.0163855061091, 0.9824051281495, 1.04908985213703, 1.01102197270052, // denom2
        1.0074050784795, 1.0018402553389, 0.94654127396428, 0.94051319172603  // denom3
#elif (NUM_CODENS == 11)
        1, 3.7812662079919, -2.0909314246215, -2.20498028081056,       // denom1
        1, -1.12941374204415, -0.326323026533379, 1.27424829228824,    // denom2
        1, -2.38787095655157, 1.29218756210116, 2.71081179958469       // denom3
#endif
};

//  3-way linear model coefficients
static const double linear_coeffs_3way_9CODENs[17] =
    {
#if (NUM_CODENS == 9 || NUM_CODENS == 10)
        -0.587100397752123,    // (Intercept)
        -7.93027271403258e-08, // XnY_ctx
        4.98481429879899e-06,  // N_diff_obj_section
        -1.63903694463017e-05, // N_mut2_ctx
        -4.93356600478016e-06, // N_diff_obj
        -14843.5350131617,     // denom1
        306478.144745683,      // denom2
        -290254.345158429,     // denom3
        -1781.00695309133,     // XnY_ctx:denom1
        -1091.99589131905,     // N_diff_obj_section:denom1
        327.6575530681,        // N_mut2_ctx:denom1
        47395.255680747,       // XnY_ctx:denom2
        47798.7112013701,      // N_diff_obj_section:denom2
        -35574.6303634354,     // N_mut2_ctx:denom2
        -45398.423188085,      // XnY_ctx:denom3
        -46461.4920783106,     // N_diff_obj_section:denom3
        35045.3583333325,      // N_mut2_ctx:denom3
#elif (NUM_CODENS == 11)
        -62.5912401461314,     // (Intercept)
        -4.63056017218561e-07, // XnY_ctx
        -5.94542251217161e-06, // N_diff_obj_section
        -1.88405508374364e-05, // N_mut2_ctx
        1.6482753400726e-05,   // N_diff_obj
        -121.674835624081,     // denom1
        -15.908243723878,      // denom2
        135.912637012686,      // denom3
        84.3989913963196,      // XnY_ctx:denom1
        16.1882290423456,      // N_diff_obj_section:denom1
        119.887931220681,      // N_mut2_ctx:denom1
        -21.4985877863572,     // XnY_ctx:denom2
        -18.6403364210245,     // N_diff_obj_section:denom2
        -12.3458956197811,     // N_mut2_ctx:denom2
        -0.305038962888618,    // XnY_ctx:denom3
        2.98700710443235,      // N_diff_obj_section:denom3
        -6.58634757454746,     // N_mut2_ctx:denom3
#endif
};
// ANI>95 subset of 3-way linear model coefficients
static const double N95linear_coeffs_3way_9CODENs[17] =
    {
#if (NUM_CODENS == 9 || NUM_CODENS == 10)
        10.7394186800999,      // (Intercept)
        -2.92835943410177e-08, // XnY_ctx
        1.66112325940088e-05,  // N_diff_obj_section
        -3.48757817098253e-05, // N_mut2_ctx
        -1.82725607919882e-05, // N_diff_obj
        3151.11534599359,      // denom1
        -1300.37794535582,     // denom2
        -1824.96575177919,     // denom3
        -21820.9227396359,     // XnY_ctx:denom1
        -30943.09511694,       // N_diff_obj_section:denom1
        10100.3465983902,      // N_mut2_ctx:denom1
        25441.8492675468,      // XnY_ctx:denom2
        35763.8363380932,      // N_diff_obj_section:denom2
        -3901.01958522226,     // N_mut2_ctx:denom2
        -3669.6748819999,      // XnY_ctx:denom3
        -4909.59256464219,     // N_diff_obj_section:denom3
        -5653.59698386914,     // N_mut2_ctx:denom3
#elif (NUM_CODENS == 11)
        -1.05834423005556,     // (Intercept)
        -4.11681306702505e-08, // XnY_ctx
        -1.83412193040086e-05, // N_diff_obj_section
        3.68710040747636e-05,  // N_mut2_ctx
        5.22422757620718e-06,  // N_diff_obj
        30.5002312286384,      // denom1
        -150.963243575499,     // denom2
        121.876576325283,      // denom3
        0.221420061509286,     // XnY_ctx:denom1
        0.563241174219071,     // N_diff_obj_section:denom1
        0.539838644833922,     // N_mut2_ctx:denom1
        1.19759583508643,      // XnY_ctx:denom2
        0.0691552803067527,    // N_diff_obj_section:denom2
        -1.28306951996752,     // N_mut2_ctx:denom2
        -0.359874242747578,    // XnY_ctx:denom3
        -0.00075252953123542,  // N_diff_obj_section:denom3
        -0.525360651962046,    // N_mut2_ctx:denom3
#endif

};
/*
//--------------------------------------------------------------//
// model 1: 3-way linear model parameters (NUM_CODENS==11 only)

// optimized parameters for 3 denominator of 3-way linear model
static const double T11opt_denom_params4x3[12] =
    {
        1.0178887, 1.0135146, 1.0302181, 1.1205122, // denom1
        1.0136223, 0.9803396, 1.0021964, 0.9839395, // denom2
        0.9988935, 0.9690589, 0.9888171, 0.9817371  // denom3
};
// ANI>95 subset parameters
static const double T11N95opt_denom_params4x3[12] =
    {
        1.0145656, 1.2708196, 1.2169829, 0.3904416,  // denom1
        0.8147527, 1.3389150, 2.0581364, -1.4624736, // denom2
        1.6877678, 2.4390067, 0.7878482, 0.6697437   // denom3
};

//  3-way linear model coefficients
static const double linear_coeffs_3way_11CODENs[17] =
    {
        6.86905236142636,      // (Intercept)
        -5.85271800517096e-08, // XnY_ctx
        7.43351941016881e-06,  // N_diff_obj_section
        -1.95859228059721e-05, // N_mut2_ctx
        -8.44336369105416e-06, // N_diff_obj
        -11726.5478433625,     // denom1
        -111109.682632437,     // denom2
        121000.895632267,      // denom3
        -741.538955661909,     // XnY_ctx:denom1
        -1894.53741972785,     // N_diff_obj_section:denom1
        2585.81964981999,      // N_mut2_ctx:denom1
        -7744.5216892858,      // XnY_ctx:denom2
        -18190.2514369096,     // N_diff_obj_section:denom2
        23306.3090989115,      // N_mut2_ctx:denom2
        8352.82760552325,      // XnY_ctx:denom3
        19782.4294837312,      // N_diff_obj_section:denom3
        -25518.2493259614,     // N_mut2_ctx:denom3
};
// ANI>95 subset of 3-way linear model coefficients
static const double N95linear_coeffs_3way_11CODENs[17] =
    {
        -1.69173585242774,     // (Intercept)
        -3.64716315423252e-08, // XnY_ctx
        2.85880264052975e-05,  // N_diff_obj_section
        -5.5597924264054e-05,  // N_mut2_ctx
        -2.96208613804271e-05, // N_diff_obj
        1106.96550282022,      // denom1
        -83.197177767225,      // denom2
        -1671.88720051458,     // denom3
        152.75590068948,       // XnY_ctx:denom1
        212.734805232395,      // N_diff_obj_section:denom1
        208.388739816682,      // N_mut2_ctx:denom1
        -0.27860125939004,     // XnY_ctx:denom2
        2.27065806912198,      // N_diff_obj_section:denom2
        -10.2977013026837,     // N_mut2_ctx:denom2
        -250.681324281127,     // XnY_ctx:denom3
        -403.751360812353,     // N_diff_obj_section:denom3
        -137.96628693833,      // N_mut2_ctx:denom3

};
*/

static const double af_ANIb_denom_params5x3[15] =
    {
        2.36775520439829, 2.17036643548818, -0.0977124271437818, -0.251063182363729, 2.19491057829433, // denom1
        1.41850238263675, 2.15224961402778, -2.18587619877498, -2.21943802445092, 0.418474287083788,   // denom2
        0.990424249133548, 0.502163745411501, 5.7787290524088, 0.626839640715682, 0.638706160147675    // denom3
};

static const double af_ANIb_3way_CODENs[21] =
    {
        0.127498260534708,                  // (Intercept)
        -1.51730137145957e-06,              // XnY_ctx
        9.48114616179122e-05,               // N_diff_obj_section
        -0.000257916039708839,              // N_mut2_ctx
        -7.85253629216324e-05,              // N_diff_obj
        1.96482915215685e-06,               // X_ctx
        -1475.58996322769,                  // denom1
        26.4698377745167,                   // denom2
        529.387641078915,                   // denom3
        42.0602258298213,                   // XnY_ctx:denom1
        92.3637001980455,                   // N_diff_obj_section:denom1
        -306.00755719743,                   // N_mut2_ctx:denom1
        -86.3482066935661,                  // N_diff_obj:denom1
        2.16630645286942,                   // XnY_ctx:denom2
        -0.674746089642152,                 // N_diff_obj_section:denom2
        48.1147600903939,                   // N_mut2_ctx:denom2
        3.39101087242059,                   // N_diff_obj:denom2
        -15.5481688297158,                  // XnY_ctx:denom3
        -26.7068442368566,                  // N_diff_obj_section:denom3
        18.5142668221724,                   // N_mut2_ctx:denom3
        20.4430437567379,                   // N_diff_obj:denom3
};

typedef struct
{
    uint32_t XnY_ctx;
    uint32_t X_ctx;
    uint32_t N_diff_obj_section;
    uint32_t N_mut2_ctx;
    uint32_t N_diff_obj;
} ani_features_t;

#define EPSILON (1e-8)

/* Universal sketch-size-normalized ANI model.
 * The installed minco calibration is trained at 10,000 bottom-k contexts.
 * Count-like inputs from other sketch sizes are normalized to that reference
	 * scale; S=10000 is exactly the trained behavior.  ani_model_compat_filter_shift is kept as
	 * a compatibility fallback for old stat files that do not carry minco metadata.
 */
enum { ANI_MODEL_REFERENCE_FILTER_SHIFT = 8 };
enum { ANI_MODEL_REFERENCE_SKETCH_SIZE = 10000 };
extern int ani_model_compat_filter_shift;
extern uint32_t ani_model_target_sketch_size;

static inline double ani_model_fold_scale(void)
{
    if (ani_model_target_sketch_size > 0)
        return (double)ANI_MODEL_REFERENCE_SKETCH_SIZE /
               (double)ani_model_target_sketch_size;
    const int delta = ani_model_compat_filter_shift - ANI_MODEL_REFERENCE_FILTER_SHIFT;
    if (delta < -60)
        return ldexp(1.0, -60);
    if (delta > 60)
        return ldexp(1.0, 60);
    return ldexp(1.0, delta);
}

static inline double lm3ways_dist_from_values_core(double XnY_ctx,
                                                   double N_diff_obj_section,
                                                   double N_mut2_ctx,
                                                   double N_diff_obj,
                                                   const double p[12],
                                                   const double coeffs[17])
{
    // compute 3 denomonators
    double denom1 = 1 / (p[0] * XnY_ctx + p[1] * N_diff_obj_section + p[2] * N_mut2_ctx + p[3] * N_diff_obj + EPSILON);
    double denom2 = 1 / (p[4] * XnY_ctx + p[5] * N_diff_obj_section + p[6] * N_mut2_ctx + p[7] * N_diff_obj + EPSILON);
    double denom3 = 1 / (p[8] * XnY_ctx + p[9] * N_diff_obj_section + p[10] * N_mut2_ctx + p[11] * N_diff_obj + EPSILON);
    // compute 3-way linear model
    double dist = coeffs[0] +
                  coeffs[1] * XnY_ctx +
                  coeffs[2] * N_diff_obj_section +
                  coeffs[3] * N_mut2_ctx +
                  coeffs[4] * N_diff_obj +
                  coeffs[5] * denom1 +
                  coeffs[6] * denom2 +
                  coeffs[7] * denom3 +
                  coeffs[8] * XnY_ctx * denom1 +
                  coeffs[9] * N_diff_obj_section * denom1 +
                  coeffs[10] * N_mut2_ctx * denom1 +
                  coeffs[11] * XnY_ctx * denom2 +
                  coeffs[12] * N_diff_obj_section * denom2 +
                  coeffs[13] * N_mut2_ctx * denom2 +
                  coeffs[14] * XnY_ctx * denom3 +
                  coeffs[15] * N_diff_obj_section * denom3 +
                  coeffs[16] * N_mut2_ctx * denom3;
    return dist;
}

static inline double lm3ways_dist_from_features_core(ani_features_t *features, const double p[12], const double coeffs[17])
{
    return lm3ways_dist_from_values_core((double)features->XnY_ctx,
                                         (double)features->N_diff_obj_section,
                                         (double)features->N_mut2_ctx,
                                         (double)features->N_diff_obj,
                                         p, coeffs);
}

// 2. 3-ways linear model for distance
static inline double lm3ways_dist_from_features(ani_features_t *features)
{
    const double scale = ani_model_fold_scale();
    const double XnY_ctx = (double)features->XnY_ctx * scale;
    const double N_diff_obj_section = (double)features->N_diff_obj_section * scale;
    const double N_mut2_ctx = (double)features->N_mut2_ctx * scale;
    const double N_diff_obj = (double)features->N_diff_obj * scale;
    if (XnY_ctx <= 0.0)
        return 1; // if no ctx, return 1
    else if (N_diff_obj <= 0.0)
        return 0; // if no diff obj, return 0
    // use optimized parameters lm prediction
    double dist = lm3ways_dist_from_values_core(XnY_ctx, N_diff_obj_section, N_mut2_ctx, N_diff_obj,
                                                opt_denom_params4x3, linear_coeffs_3way_9CODENs);
    if (dist < 0.05)
        dist = lm3ways_dist_from_values_core(XnY_ctx, N_diff_obj_section, N_mut2_ctx, N_diff_obj,
                                             N95opt_denom_params4x3, N95linear_coeffs_3way_9CODENs);
    if (dist < 0)
        dist = 0;
    return dist;
}
// 3. 3-ways linear model for af_ANIb
static inline double lm3ways_af_ANIb_from_values_core(double XnY_ctx,
                                                      double N_diff_obj_section,
                                                      double N_mut2_ctx,
                                                      double N_diff_obj,
                                                      double X_ctx,
                                                      const double p[15],
                                                      const double coeffs[21])
{
    // compute 3 denomonators
    double denom1 = 1 / (p[0] * XnY_ctx + p[1] * N_diff_obj_section + p[2] * N_mut2_ctx + p[3] * N_diff_obj + p[4] * X_ctx + EPSILON);
    double denom2 = 1 / (p[5] * XnY_ctx + p[6] * N_diff_obj_section + p[7] * N_mut2_ctx + p[8] * N_diff_obj + p[9] * X_ctx + EPSILON);
    double denom3 = 1 / (p[10] * XnY_ctx + p[11] * N_diff_obj_section + p[12] * N_mut2_ctx + p[13] * N_diff_obj + p[14] * X_ctx + EPSILON);
    
    // compute 3-way linear model
    double af_ANIb = coeffs[0] +
                  coeffs[1] * XnY_ctx +
                  coeffs[2] * N_diff_obj_section +
                  coeffs[3] * N_mut2_ctx +
                  coeffs[4] * N_diff_obj +
                  coeffs[5] * X_ctx +
                  coeffs[6] * denom1 +
                  coeffs[7] * denom2 +
                  coeffs[8] * denom3 +
                  coeffs[9] * XnY_ctx * denom1 +
                  coeffs[10] * N_diff_obj_section * denom1 +
                  coeffs[11] * N_mut2_ctx * denom1 +
                  coeffs[12] * N_diff_obj * denom1 +
                  coeffs[13] * XnY_ctx * denom2 +
                  coeffs[14] * N_diff_obj_section * denom2 +
                  coeffs[15] * N_mut2_ctx * denom2 +
                  coeffs[16] * N_diff_obj * denom2 +
                  coeffs[17] * XnY_ctx * denom3 +
                  coeffs[18] * N_diff_obj_section * denom3 +
                  coeffs[19] * N_mut2_ctx * denom3 +
                  coeffs[20] * N_diff_obj * denom3;

    return af_ANIb;
}

static inline double lm3ways_af_ANIb_from_features_core(ani_features_t *features, const double p[15], const double coeffs[21])
{
    return lm3ways_af_ANIb_from_values_core((double)features->XnY_ctx,
                                            (double)features->N_diff_obj_section,
                                            (double)features->N_mut2_ctx,
                                            (double)features->N_diff_obj,
                                            (double)features->X_ctx,
                                            p, coeffs);
}

static inline double lm3ways_af_ANIb_from_features(ani_features_t *features)
{
    // use optimized parameters lm prediction
    const double scale = ani_model_fold_scale();
    double af_ANIb = lm3ways_af_ANIb_from_values_core((double)features->XnY_ctx * scale,
                                                      (double)features->N_diff_obj_section * scale,
                                                      (double)features->N_mut2_ctx * scale,
                                                      (double)features->N_diff_obj * scale,
                                                      (double)features->X_ctx * scale,
                                                      af_ANIb_denom_params5x3, af_ANIb_3way_CODENs);
    return af_ANIb;
}   




#endif
