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