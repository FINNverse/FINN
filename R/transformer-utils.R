## Transformer modules
library(torch)
library(R6)
library(glue)


canonical_mask = function(mask = NULL, mask_name = NULL, other_type = NULL, other_name = NULL, target_type = NULL, check_other = TRUE) {
  if(!is.null(mask)) {
    mask_dtype = mask$dtype
    mask_is_float = torch_is_floating_point(mask)

    if (mask_dtype != torch_bool() && (!mask_is_float)) {
      stop(glue::glue("only bool and floating types of {mask_name} are supported"))
    }

    if(!is.null(check_other) && !is.null(other_type)) {
      if(mask_dtype != other_type) {
        warning(glue::glue("Support for mismatched {mask_name} and {other_name} is deprecated. Use same type for both instead. "))
      }
    }

    if(!mask_is_float) {
      mask = torch_zeros_like(mask, dtype=target_type)$masked_fill_(mask, -Inf)
    }
  }
  return(mask)
}

generate_square_subsequent_mask = function(sz, device = NULL, dtype = NULL) {
  return( torch_triu(torch_full(c(sz, sz), -Inf, dtype=dtype, device=device), diagonal = 1L))
}

detect_is_causal_mask = function(mask, is_causal = NULL, size = NULL) {

  make_causal = FALSE

  if(!is.null(is_causal)) {
    if(is_causal) {
      make_causal = TRUE
    } else {
      make_causal = FALSE
    }
  }

  if(is.null(is_causal) && !is.null(mask)) {
    sz = if(!is.null(size)) size else mask$size(-2)
    causal_comparison = generate_square_subsequent_mask(sz, device=mask$device, dtype=mask$dtype)

    if(all(mask$size() & causal_comparison$size())) {
      make_causal = (mask == causal_comparison)$all() |> as.numeric() |> as.logical()
    } else
      make_causal = FALSE
  }
  return(make_causal)
}

TransformerEncoderLayer = nn_module(
  classname = "TransformerEncoderLayer",
  initialize = function(d_model,
                        nhead,
                        dim_feedforward = 2048,
                        dropout = 0.1,
                        activation = torch::nnf_relu,
                        layer_norm_eps = 1e-5,
                        batch_first = FALSE,
                        norm_first = FALSE,
                        bias = TRUE,
                        device="cpu",
                        dtype=torch::torch_float32()) {

    self$self_attn = torch::nn_multihead_attention(
      embed_dim = d_model,
      num_heads = nhead,
      dropout = dropout,
      bias = bias,
      batch_first = batch_first
    )

    self$linear1 = nn_linear(d_model, dim_feedforward, bias=bias)
    self$dropout = nn_dropout(dropout)
    self$linear2 = nn_linear(dim_feedforward, d_model, bias=bias)
    self$norm_first = norm_first
    self$norm1 = nn_layer_norm(d_model, eps=layer_norm_eps)
    self$norm2 = nn_layer_norm(d_model, eps=layer_norm_eps)
    self$dropout1 = nn_dropout(dropout)
    self$dropout2 = nn_dropout(dropout)
    self$activation = activation
    self$save_attn_weights = FALSE
    self$attn_weights = NA
  },

  forward = function(src,
                     src_mask = NULL,
                     src_key_padding_mask = NULL,
                     is_causal = FALSE) {

    src_key_padding_mask = canonical_mask(mask = src_key_padding_mask,
                                          mask_name = "src_key_padding_mask",
                                          other_type = NULL,
                                          other_name = "",
                                          #target_type = src$dtype,
                                          check_other = FALSE)

    src_mask = canonical_mask(mask = src_mask,
                              mask_name = "src_mask",
                              other_type = NULL,
                              other_name = "",
                              target_type = src$dtype,
                              check_other = FALSE)

    x = src
    if (self$norm_first) {
      x = x + self$sa_block(
        x = self$norm1(x), att_mask = src_mask, key_padding_mask = src_key_padding_mask, is_causal = is_causal)
      x = x + self$ff_block(self$norm2(x))
    } else {
      x = self$norm1(
        x + self$sa_block(x = x, att_mask = src_mask, key_padding_mask = src_key_padding_mask, is_causal = is_causal)
      )
      x = self$norm2(x + self$ff_block(x))
    }
    return(x)
  },

  sa_block = function(x,
                      att_mask = NULL,
                      key_padding_mask = NULL,
                      is_causal = FALSE
  ) {
    x = self$self_attn(
      x,
      x,
      x,
      attn_mask = att_mask,
      key_padding_mask = key_padding_mask,
      need_weights = TRUE
    )
    #if(self$save_attn_weights) {
    self$attn_weights = x[[2]]
    #}
    return( self$dropout1(x[[1]]) )
  },

  ff_block = function(x) {
    x = self$linear2(self$dropout(self$activation(self$linear1(x))))
    return(self$dropout2(x))
  }
)

get_seq_len = function(src, batch_first) {
  src_size = src$size()
  if(length(src_size) == 2) return(src_size[1])
  else {
    seq_len_pos = if(batch_first) 2 else 1
    return(seq_len_pos)
  }
}

TransformerEncoder = nn_module(
  classname = "Transformer",
  initialize = function(encoder_layer = TransformerEncoderLayer(d_model = 512, nhead = 4L),
                        num_layers,
                        norm = NULL,
                        mask_check = TRUE) {

    for(i in 1:num_layers) {
      self[[paste0("layers_", i)]] = encoder_layer$clone(deep = TRUE)
    }
    self$num_layers = num_layers
    self$norm = norm
    self$mask_check = mask_check
  },

  forward = function(src,
                     mask = NULL,
                     src_key_padding_mask = NULL,
                     is_causal = NULL) {

    src_key_padding_mask = canonical_mask(mask = src_key_padding_mask,
                                          mask_name = "src_key_padding_mask",
                                          other_type = NULL,
                                          other_name = "",
                                          #target_type = src$dtype,
                                          check_other = FALSE)

    src_mask = canonical_mask(mask = mask,
                              mask_name = "src_mask",
                              other_type = NULL,
                              other_name = "",
                              target_type = src$dtype,
                              check_other = FALSE)

    output = src
    first_layer = self$layers_1
    src_key_padding_mask_for_layers = src_key_padding_mask
    batch_first = first_layer$self_attn$batch_first

    seq_len = get_seq_len(src, batch_first)
    is_causal = detect_is_causal_mask(mask, is_causal)

    for(i in 1:self$num_layers) {
      output = self[[paste0("layers_", i)]](
        output,
        src_mask = mask,
        is_causal = is_causal,
        src_key_padding_mask=src_key_padding_mask_for_layers
      )
    }
    if (!is.null(self$norm )) {
      output = self$norm(output)
    }
    return(output)
  }
)

nbinom_torch = function(pred, true, theta) {
  eps = 0.0001
  pred = torch::torch_exp(pred)
  theta_tmp = 1.0/(torch::nnf_softplus(theta)+eps)
  probs = torch::torch_clamp(1.0 - theta_tmp/(theta_tmp+pred)+eps, 0.0, 1.0-eps)
  total_count = theta_tmp
  value = true
  logits = torch::torch_log(probs) - torch::torch_log1p(-probs)
  log_unnormalized_prob <- total_count * torch::torch_log(torch::torch_sigmoid(-logits)) + value * torch::torch_log(torch::torch_sigmoid(logits))
  log_normalization <- -torch::torch_lgamma(total_count + value) + torch::torch_lgamma(1 + value) + torch::torch_lgamma(total_count)
  log_normalization <- torch::torch_where(total_count + value == 0, torch::torch_tensor(0, dtype = log_normalization$dtype, device = pred$device), log_normalization)
  return( - (log_unnormalized_prob - log_normalization))
}


PositionalEncoding = nn_module("PositionalEncoding",
                               initialize = function(emb_dim = 50L, dropout=0.1, max_len=50L) {
                                 super$initialize()
                                 self$dropout = nn_dropout(p = dropout)
                                 pe = torch_zeros(max_len, emb_dim)
                                 position = torch_arange(1, max_len, dtype=torch_float32())$unsqueeze(2L)
                                 div_term = torch_exp(torch_arange(1, emb_dim, 2)$float() * (-log(10000.0) / emb_dim))

                                 ind = which((1:emb_dim) %% 2 == 0, arr.ind = TRUE)
                                 pe[, ind] =  if(length(pe[, ind]$shape) == 1) torch_sin(position * div_term)[,1] else torch_sin(position * div_term)
                                 ind = which((1:emb_dim) %% 2 != 0, arr.ind = TRUE)
                                 pe[, ind] = if(length(pe[, ind]$shape) == 1) torch_cos(position * div_term)[,1] else torch_cos(position * div_term)

                                 pe = pe$unsqueeze(1L)
                                 self$register_buffer("pe", pe)
                               },

                               forward = function(x) {
                                 x = x + self$pe[,1:x$shape[2],]
                                 return(self$dropout(x))
                               }
)
