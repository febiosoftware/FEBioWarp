#pragma once

class Image;

void blur_image_2d(Image& trg, Image& src, float d);
void blur_image(Image& trg, Image& src, float d);
