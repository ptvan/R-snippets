library(keras)

model <- application_vgg16(weights = 'imagenet', include_top = FALSE)

# this call requires the Python Imaging Library fork (`pip3 install pillow`)
img <- image_load("Flicker8k_Dataset/", target_size = c(224,224))

x <- image_to_array(img)
x <- array_reshape(x, c(1, dim(x)))
x <- imagenet_preprocess_input(x)

features <- model %>% predict