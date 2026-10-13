#' @importFrom stats coef cor.test dfbetas
NULL

#' @title Getting DFBETA values
#'
#' @description
#' Returns the raw, signed, and signed and ranked DFBETA values calculated from
#' the given linear regression models representing a mediation model.
#'
#' @param model_a Regression model of X predicting M
#' @param model_bc Regression model of X and M predicting Y
#' @param X (character) X variable in model
#' @param M (character) M variable in model
#'
#' @return Named list of DFBETA values
#' @export
get_dfbetas <- function(model_a, model_bc, X, M) {
    dfb_a_matrix <- stats::dfbetas(model_a)
    dfb_bc_matrix <- stats::dfbetas(model_bc)

    dfb_a <- dfb_a_matrix[, X]
    dfb_b <- dfb_bc_matrix[, M]

    dfb_a_signed <- dfb_a * sign(stats::coef(model_a)[X])
    dfb_b_signed <- dfb_b * sign(stats::coef(model_bc)[M])

    dfb_a_rank_signed <- rank(dfb_a_signed) / length(dfb_a_signed)
    dfb_b_rank_signed <- rank(dfb_b_signed) / length(dfb_b_signed)

    return(list(
        "dfb_a" = dfb_a
        , "dfb_b" = dfb_b
        , "dfb_a_path_aligned" = dfb_a_signed
        , "dfb_b_path_aligned" = dfb_b_signed
        , "dfb_a_ranked_path_aligned" = dfb_a_rank_signed
        , "dfb_b_ranked_path_aligned" = dfb_b_rank_signed
    ))
}

#' @title Getting cross-path DFBETA correlations
#'
#' @description
#' Returns cross-path DFBETA correlations calculated from the given DFBETA values.
#'
#' @param dfbetas Named list of values returned from function `get_dfbetas`
#'
#' @returns Named list of DFBETA cross-path correlations and p values
#' @export
get_dfbetas_corrs <- function(dfbetas) {
    dfb_corr_raw <- stats::cor.test(
        dfbetas$dfb_a
        , dfbetas$dfb_b
    )
    dfb_corr_raw_signed <- stats::cor.test(
        dfbetas$dfb_a_path_aligned
        , dfbetas$dfb_b_path_aligned
    )
    dfb_corr_rank <- stats::cor.test(
        dfbetas$dfb_a
        , dfbetas$dfb_b
        , method = "spearman"
    )
    dfb_corr_rank_signed <- stats::cor.test(
        dfbetas$dfb_a_path_aligned
        , dfbetas$dfb_b_path_aligned
        , method = "spearman"
    )

    return(list(
        "dfbeta_corr" = dfb_corr_raw$estimate
        , "dfbeta_path_aligned_corr" = dfb_corr_raw_signed$estimate
        , "p_value" = dfb_corr_raw$p.value
        , "dfbeta_ranked_corr" = dfb_corr_rank$estimate
        , "dfbeta_ranked_path_aligned_corr" = dfb_corr_rank_signed$estimate
        , "p_value_ranked_dfbeta" = dfb_corr_rank$p.value
    ))
}

#' @title Getting IIM correlations
#'
#' @description
#' Calculates the individual-index of mediation (IIM) from the given DFBETA values
#' and correlates between the IIM and given moderator data.
#'
#' @param data dataframe
#' @param dfbetas named list of dfbetas from `get_dfbetas`
#' @param W_cols vector of column names of moderator values to calculate iims for
#'
#' @returns Dataframe where each row contain IIM correlations with a single moderator
#' @export
get_iim_corrs <- function(data, dfbetas, W_cols) {
    iim_raw_add <- dfbetas$dfb_a + dfbetas$dfb_b
    iim_raw_minus <- dfbetas$dfb_a - dfbetas$dfb_b
    iim_signed_add <- dfbetas$dfb_a_signed + dfbetas$dfb_b_signed
    iim_signed_minus <- dfbetas$dfb_a_signed - dfbetas$dfb_b_signed
    iim_rank_add <- dfbetas$dfb_a_rank_signed + dfbetas$dfb_b_rank_signed
    iim_rank_minus <- dfbetas$dfb_a_rank_signed - dfbetas$dfb_b_rank_signed

    result <- list()

    for (col in W_cols) {
        iim_corr_raw_add <- stats::cor.test(iim_raw_add, data[[col]])
        iim_corr_raw_minus <- stats::cor.test(iim_raw_minus, data[[col]])
        iim_corr_signed_add <- stats::cor.test(iim_signed_add, data[[col]])
        iim_corr_signed_minus <- stats::cor.test(iim_signed_minus, data[[col]])
        iim_corr_rank_add <- stats::cor.test(iim_rank_add, data[[col]])
        iim_corr_rank_minus <- stats::cor.test(iim_rank_minus, data[[col]])

        dfb_a_signed_corr <- stats::cor.test(dfbetas$dfb_a_signed, data[[col]])
        dfb_b_signed_corr <- stats::cor.test(dfbetas$dfb_b_signed, data[[col]])

        result[[col]] <- data.frame(
            "moderator" = col
            , "iim_add" = iim_corr_raw_add$estimate
            , "p_val_add" = iim_corr_raw_add$p.value
            , "iim_minus" = iim_corr_raw_minus$estimate
            , "p_val_minus" = iim_corr_raw_minus$p.value
            , "iim_add_sign_adjusted_dfbetas" = iim_corr_signed_add$estimate
            , "p_val_add_path_aligned_dfbetas" = iim_corr_signed_add$p.value
            , "iim_minus_path_aligned_dfbetas" = iim_corr_signed_minus$estimate
            , "p_val_minus_path_aligned_dfbetas" = iim_corr_signed_minus$p.value
            , "iim_add_ranked_path_aligned_dfbetas" = iim_corr_rank_add$estimate
            , "p_val_iim_add_ranked_path_aligned_dfbetas" = iim_corr_rank_add$p.value
            , "iim_minus_ranked_path_aligned_dfbetas" = iim_corr_rank_minus$estimate
            , "p_val_iim_minus_ranked_path_aligned_dfbetas" = iim_corr_rank_minus$p.value
            , "dfb_a_path_aligned_corr" = dfb_a_signed_corr$estimate
            , "p_val_dfb_a_path_aligned_corr" = dfb_a_signed_corr$p.value
            , "dfb_b_path_aligned_corr" = dfb_b_signed_corr$estimate
            , "p_val_dfb_b_path_aligned_corr" = dfb_b_signed_corr$p.value
        )
    }

    result <- as.data.frame(do.call(rbind, result))
    return(result)
}
