library(keras)

### MNIST classification
# load data
mnist <- dataset_mnist()

x_train <- mnist$train$x
y_train <- mnist$train$y
x_test <- mnist$test$x
y_test <- mnist$test$y

# reshape
x_train <- array_reshape(x_train, c(nrow(x_train), 784))
x_test <- array_reshape(x_test, c(nrow(x_test), 784))

# rescale
x_train <- x_train / 255
x_test <- x_test / 255

# one-hot encode
y_train <- to_categorical(y_train, 10)
y_test <- to_categorical(y_test, 10)

model <- keras_model_sequential() 

model %>% 
  layer_dense(units = 256, activation = 'relu', input_shape = c(784)) %>% 
  layer_dropout(rate = 0.4) %>% 
  layer_dense(units = 128, activation = 'relu') %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 10, activation = 'softmax')

summary(model)

model %>% compile(
  loss = 'categorical_crossentropy',
  optimizer = optimizer_rmsprop(),
  metrics = c('accuracy')
)

history <- model %>% fit(
  x_train, y_train, 
  epochs = 30, batch_size = 128, 
  validation_split = 0.2
)

model %>% evaluate(x_test, y_test)

model %>% predict_classes(x_test)

### VGG16
model <- application_vgg16(weights = 'imagenet', include_top = FALSE)

# this call requires the Python Imaging Library fork (`pip3 install pillow`)
img <- image_load("Flicker8k_Dataset/", target_size = c(224,224))

x <- image_to_array(img)
x <- array_reshape(x, c(1, dim(x)))
x <- imagenet_preprocess_input(x)

features <- model %>% predict

### Generative Adversarial Network

# generator custom model
generator <-
  function(name = NULL) {
    keras_model_custom(name = name, function(self) {
      
      self$fc1 <- layer_dense(units = 7 * 7 * 64, use_bias = FALSE)
      self$batchnorm1 <- layer_batch_normalization()
      self$leaky_relu1 <- layer_activation_leaky_relu()
      self$conv1 <-
        layer_conv_2d_transpose(
          filters = 64,
          kernel_size = c(5, 5),
          strides = c(1, 1),
          padding = "same",
          use_bias = FALSE
        )
      self$batchnorm2 <- layer_batch_normalization()
      self$leaky_relu2 <- layer_activation_leaky_relu()
      self$conv2 <-
        layer_conv_2d_transpose(
          filters = 32,
          kernel_size = c(5, 5),
          strides = c(2, 2),
          padding = "same",
          use_bias = FALSE
        )
      self$batchnorm3 <- layer_batch_normalization()
      self$leaky_relu3 <- layer_activation_leaky_relu()
      self$conv3 <-
        layer_conv_2d_transpose(
          filters = 1,
          kernel_size = c(5, 5),
          strides = c(2, 2),
          padding = "same",
          use_bias = FALSE,
          activation = "tanh"
        )
      
      function(inputs, mask = NULL, training = TRUE) {
        self$fc1(inputs) %>%
          self$batchnorm1(training = training) %>%
          self$leaky_relu1() %>%
          k_reshape(shape = c(-1, 7, 7, 64)) %>%
          self$conv1() %>%
          self$batchnorm2(training = training) %>%
          self$leaky_relu2() %>%
          self$conv2() %>%
          self$batchnorm3(training = training) %>%
          self$leaky_relu3() %>%
          self$conv3()
      }
    })
  }

# discriminator custom model
discriminator <-
  function(name = NULL) {
    keras_model_custom(name = name, function(self) {
      
      self$conv1 <- layer_conv_2d(
        filters = 64,
        kernel_size = c(5, 5),
        strides = c(2, 2),
        padding = "same"
      )
      self$leaky_relu1 <- layer_activation_leaky_relu()
      self$dropout <- layer_dropout(rate = 0.3)
      self$conv2 <-
        layer_conv_2d(
          filters = 128,
          kernel_size = c(5, 5),
          strides = c(2, 2),
          padding = "same"
        )
      self$leaky_relu2 <- layer_activation_leaky_relu()
      self$flatten <- layer_flatten()
      self$fc1 <- layer_dense(units = 1)
      
      function(inputs, mask = NULL, training = TRUE) {
        inputs %>% self$conv1() %>%
          self$leaky_relu1() %>%
          self$dropout(training = training) %>%
          self$conv2() %>%
          self$leaky_relu2() %>%
          self$flatten() %>%
          self$fc1()
      }
    })
  }

# model creation
generator <- generator()
discriminator <- discriminator()

generator$call = tf$contrib$eager$defun(generator$call)
discriminator$call = tf$contrib$eager$defun(discriminator$call)

discriminator_loss <- function(real_output, generated_output) {
  real_loss <- tf$losses$sigmoid_cross_entropy(
    multi_class_labels = k_ones_like(real_output),
    logits = real_output)
  generated_loss <- tf$losses$sigmoid_cross_entropy(
    multi_class_labels = k_zeros_like(generated_output),
    logits = generated_output)
  real_loss + generated_loss
}

generator_loss <- function(generated_output) {
  tf$losses$sigmoid_cross_entropy(
    tf$ones_like(generated_output),
    generated_output)
}

discriminator_optimizer <- tf$train$AdamOptimizer(1e-4)
generator_optimizer <- tf$train$AdamOptimizer(1e-4)

# training loop
train <- function(dataset, epochs, noise_dim) {
  for (epoch in seq_len(num_epochs)) {
    start <- Sys.time()
    total_loss_gen <- 0
    total_loss_disc <- 0
    iter <- make_iterator_one_shot(train_dataset)
    
    until_out_of_range({
      batch <- iterator_get_next(iter)
      noise <- k_random_normal(c(batch_size, noise_dim))
      with(tf$GradientTape() %as% gen_tape, { with(tf$GradientTape() %as% disc_tape, {
        generated_images <- generator(noise)
        disc_real_output <- discriminator(batch, training = TRUE)
        disc_generated_output <-
          discriminator(generated_images, training = TRUE)
        gen_loss <- generator_loss(disc_generated_output)
        disc_loss <-
          discriminator_loss(disc_real_output, disc_generated_output)
      }) })
      
      gradients_of_generator <-
        gen_tape$gradient(gen_loss, generator$variables)
      gradients_of_discriminator <-
        disc_tape$gradient(disc_loss, discriminator$variables)
      
      generator_optimizer$apply_gradients(purrr::transpose(
        list(gradients_of_generator, generator$variables)
      ))
      discriminator_optimizer$apply_gradients(purrr::transpose(
        list(gradients_of_discriminator, discriminator$variables)
      ))
      
      total_loss_gen <- total_loss_gen + gen_loss
      total_loss_disc <- total_loss_disc + disc_loss
      
    })
    
    cat("Time for epoch ", epoch, ": ", Sys.time() - start, "\n")
    cat("Generator loss: ", total_loss_gen$numpy() / batches_per_epoch, "\n")
    cat("Discriminator loss: ", total_loss_disc$numpy() / batches_per_epoch, "\n\n")
    if (epoch %% 10 == 0)
      generate_and_save_images(generator,
                               epoch,
                               random_vector_for_generation)
    
  }
}

# function to save output images
generate_and_save_images <- function(model, epoch, test_input) {
  predictions <- model(test_input, training = FALSE)
  png(paste0("images_epoch_", epoch, ".png"))
  par(mfcol = c(5, 5))
  par(mar = c(0.5, 0.5, 0.5, 0.5),
      xaxs = 'i',
      yaxs = 'i')
  for (i in 1:25) {
    img <- predictions[i, , , 1]
    img <- t(apply(img, 2, rev))
    image(
      1:28,
      1:28,
      img * 127.5 + 127.5,
      col = gray((0:255) / 255),
      xaxt = 'n',
      yaxt = 'n'
    )
  }
  dev.off()
}

# actually run the generation
num_epochs <- 150
train(train_dataset, num_epochs, noise_dim)