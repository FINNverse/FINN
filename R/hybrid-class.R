########## Hyb ########
hybrid_transformer = nn_module("hybrid_transformer",
                            initialize = function(num_species = 5L,
                                                  num_env_vars =7L,
                                                  dgtl_embedder_dim = 4L,
                                                  max_len = 1000L,
                                                  emb_dim=20L,
                                                  num_heads=1L,
                                                  num_layers=1L,
                                                  dropout=0.1,
                                                  dim_feedforward = 512L
                            ) {
                              self$num_species = num_species
                              self$species_embedder = nn_embedding(num_embeddings=num_species, embedding_dim=emb_dim/2)
                              self$dgtl_embedder = nn_sequential(nn_linear(dgtl_embedder_dim, 50L),
                                                                 nn_leaky_relu(),
                                                                 nn_dropout(dropout),
                                                                 nn_linear(50L, emb_dim/2))
                              self$env_embedder = nn_sequential(nn_linear(num_env_vars, 50L),
                                                                nn_leaky_relu(),
                                                                nn_dropout(dropout),
                                                                nn_linear(50L, emb_dim))

                              self$positional_encoding = PositionalEncoding(emb_dim, dropout, max_len = max_len)
                              encoder_layer = TransformerEncoderLayer(d_model=emb_dim, nhead=num_heads,batch_first = TRUE, dim_feedforward = dim_feedforward)
                              self$transformer_encoder = TransformerEncoder(encoder_layer, num_layers=num_layers)
                              self$output = nn_sequential(nn_linear(emb_dim, 100L),
                                                          nn_leaky_relu(),
                                                          nn_dropout(dropout),
                                                          nn_linear(100, 1L))

                            },
                            forward = function(dbh, growth = NULL, trees, light, species, env) {
                              # cast patches*batch / cast at the end back
                              orig_shape1 = dbh$shape[1:2]
                              orig_shape2 = dbh$shape[3]

                              #browser()

                              if(!is.null(growth)) {
                                cohort_context = torch_cat(list(
                                  torch::torch_log(dbh$view(c(-1L, orig_shape2))$unsqueeze(3L)+1.0),
                                  growth$view(c(-1L, orig_shape2))$unsqueeze(3L),
                                  torch::torch_log(trees$view(c(-1L, orig_shape2))$unsqueeze(3L)+1.0),
                                  light$view(c(-1L, orig_shape2))$unsqueeze(3L)
                                ), dim = 3L)
                              } else {
                                cohort_context = torch_cat(list(
                                  torch::torch_log(dbh$view(c(-1L, orig_shape2))$unsqueeze(3L)+1.0),
                                  torch::torch_log(trees$view(c(-1L, orig_shape2))$unsqueeze(3L)+1.0),
                                  light$view(c(-1L, orig_shape2))$unsqueeze(3L)
                                ), dim = 3L)
                              }
                              cohort_context_em = self$dgtl_embedder(cohort_context)
                              species_em = self$species_embedder(species$view(c(-1L, orig_shape2))$to(dtype=torch::torch_long()))
                              env_em = self$env_embedder( env$unsqueeze(2L)$`repeat`(c(1L, orig_shape1[2],1L)) )
                              env_em = env_em$view(c(-1L, species_em$shape[3]*2))
                              species_em_cont = torch::torch_cat(list(
                                species_em,
                                cohort_context_em
                              ), 3L)
                              embeddings =
                                torch::torch_cat(list(
                                  env_em$unsqueeze(2L),
                                  species_em_cont
                                ), 2L)
                              embeddings = self$positional_encoding(embeddings)

                              src_key_padding_mask = trees$view(c(-1L, orig_shape2))$lt(0.5)
                              env_mask = torch_zeros(embeddings$shape[1], 1, dtype=torch_bool(), device=embeddings$device)
                              key_padding_mask = torch_cat(list(env_mask, src_key_padding_mask), dim=2)

                              #src_mask = matrix(0L, embeddings$shape[2], embeddings$shape[2])
                              #diag(src_mask)[2:nrow(src_mask)] = 1
                              #src_mask = torch_tensor(src_mask, dtype=torch_bool(), device=self$device)

                              transformer_output = self$transformer_encoder(embeddings, src_key_padding_mask=key_padding_mask)
                              pred = self$output(transformer_output[,2:(orig_shape2+1),])$squeeze(3L)
                              return(pred$view(c(orig_shape1, orig_shape2)))
                            }
)


hybrid_DNN = nn_module("hybrid_DNN",
                               initialize = function(num_species = 5L,
                                                     num_env_vars =7L,
                                                     emb_dim=20L,
                                                     dropout=0.1,
                                                     hidden = c(50L, 50L),
                                                     ...
                               ) {
                                 self$species_embedder = nn_embedding(num_species, emb_dim)

                                 layers = list()
                                 layers[[1]] = nn_linear(num_env_vars+3L+emb_dim, hidden[1])
                                 counter = 2L
                                 for(i in 2:length(hidden)) {
                                   layers[[counter]] = nn_selu()
                                   counter = counter + 1
                                   layers[[counter]] = nn_dropout(dropout)
                                   counter = counter + 1
                                   layers[[counter]] = nn_linear(hidden[i-1], hidden[i])
                                   counter = counter + 1
                                 }
                                 layers[[counter]] = nn_selu()
                                 layers[[counter+1]] = nn_dropout(dropout)
                                 layers[[counter+2]] = nn_linear(hidden[length(hidden)], 1L)
                                 self$nn = do.call(nn_sequential, layers)

                               },
                       forward = function(dbh, growth = NULL, trees, light, species, env) {
                         if(is.null(growth)) {
                           input =
                             torch::torch_cat(
                               list(dbh$unsqueeze(4L)/100., # rescale dbh wird das noch benoetigt?
                                    (trees$unsqueeze(4L)+1.0)$log(),
                                    env$unsqueeze(2L)$unsqueeze(2L)$`repeat`(c(1L, dbh$shape[2], dbh$shape[3], 1L)),
                                    light$unsqueeze(4L)),
                               dim = 4L)
                         } else {
                           input =
                             torch::torch_cat(
                               list(dbh$unsqueeze(4L)/200., # rescale dbh wird das noch benoetigt?
                                    growth$unsqueeze(4L),
                                    (trees$unsqueeze(4L)+1.0)$log(),
                                    env$unsqueeze(2L)$unsqueeze(2L)$`repeat`(c(1L, dbh$shape[2], dbh$shape[3], 1L)),
                                    light$unsqueeze(4L)),
                               dim = 4L)
                         }

                         species_embeddings = self$species_embedder(species) # [sites, patches, cohorts, species_embeddings]
                         input_together = torch::torch_cat(list(input, species_embeddings), dim = 4L)
                         predictions = self$nn(input_together) # [sites, patches, cohorts, 1]
                         return(predictions[,,,1])
                       }

)
