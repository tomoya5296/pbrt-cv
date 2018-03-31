#include "cv_film.h"
#include <memory>

#include "imageio.h"
#include "paramset.h"


namespace pbrt {

CvFilm::CvFilm(const Point2i &resolution, const Bounds2f &cropWindow,
               std::unique_ptr<Filter> filter, Float diagonal,
               const std::string &filename, Float scale,
               Float maxSampleLuminance)
    : Film(resolution, cropWindow, std::move(filter),
           diagonal, filename, scale, maxSampleLuminance) {
    cvPixels = std::unique_ptr<CvDualPixel[]>(new CvDualPixel[croppedPixelBounds.Area()]);
}

CvDualPixel &CvFilm::GetCvPixel(const Point2i &p) {
    CHECK(InsideExclusive(p, croppedPixelBounds));
    int width = croppedPixelBounds.pMax.x - croppedPixelBounds.pMin.x;
    int offset = (p.x - croppedPixelBounds.pMin.x) +
                 (p.y - croppedPixelBounds.pMin.y) * width;
    return cvPixels[offset];
}

std::unique_ptr<CvFilmTile> CvFilm::GetCvFilmTile(const Bounds2i &sampleBounds) {
    // Bound image pixels that samples in _sampleBounds_ contribute to
    Vector2f halfPixel = Vector2f(0.5f, 0.5f);
    Bounds2f floatBounds = (Bounds2f)sampleBounds;
    Point2i p0 = (Point2i)Ceil(floatBounds.pMin - halfPixel - filter->radius);
    Point2i p1 = (Point2i)Floor(floatBounds.pMax - halfPixel + filter->radius) +
                 Point2i(1, 1);
    Bounds2i tilePixelBounds = Intersect(Bounds2i(p0, p1), croppedPixelBounds);
    return std::unique_ptr<CvFilmTile>(new CvFilmTile(
        tilePixelBounds, filter->radius, filterTable, filterTableWidth,
        maxSampleLuminance));
}

void CvFilm::MergeFilmTile(std::unique_ptr<CvFilmTile> tile) {
    ProfilePhase p(Prof::MergeFilmTile);
    VLOG(1) << "Merging film tile " << tile->GetPixelBounds();
    std::lock_guard<std::mutex> lock(mutex);
    for (Point2i pixel : tile->GetPixelBounds()) {
        // Merge _pixel_ into _Film::pixels_
        const CvDualPixel &tilePixel = tile->GetPixel(pixel);
        CvDualPixel &mergePixel = GetCvPixel(pixel);
        mergePixel.AddPixel(tilePixel);
        mergePixel.filterWeightSum += tilePixel.filterWeightSum;
    }
}

void CvFilm::WriteImage(Float splatScale, int samplesPerPixel) {
    
    // Convert image to RGB and compute final pixel values
    LOG(INFO) <<
        "Converting image to RGB and computing final weighted pixel values";
    std::unique_ptr<Float[]> rgb1(new Float[3 * croppedPixelBounds.Area()]);
    std::unique_ptr<Float[]> rgb1Squared(new Float[3 * croppedPixelBounds.Area()]);
    std::unique_ptr<Float[]> rgb2(new Float[3 * croppedPixelBounds.Area()]);
    std::unique_ptr<Float[]> rgb2Squared(new Float[3 * croppedPixelBounds.Area()]);
    std::unique_ptr<Float[]> diff(new Float[3 * croppedPixelBounds.Area()]);
    std::unique_ptr<Float[]> diffSquared(new Float[3 * croppedPixelBounds.Area()]);
    std::unique_ptr<Float[]> recipPdfs(new Float[3 * croppedPixelBounds.Area()]);

    int offset = 0;
    Float avgRecipPdf = Float(0.f);
    for (Point2i p : croppedPixelBounds) {
        // Convert pixel XYZ color to RGB
        CvDualPixel &pixel = GetCvPixel(p);

        // Normalize pixel with weight sum
        Float filterWeightSum = pixel.filterWeightSum;
        Float invWt = (Float)1 / (filterWeightSum + Float(1.0e-8f));

        // Add splat value at pixel
        rgb1[3 * offset] += splatScale * pixel.L1[0] * invWt;
        rgb1[3 * offset + 1] += splatScale * pixel.L1[1] * invWt;
        rgb1[3 * offset + 2] += splatScale * pixel.L1[2] * invWt;
        rgb1Squared[3 * offset] += splatScale * pixel.L1square[0] * invWt;
        rgb1Squared[3 * offset + 1] += splatScale * pixel.L1square[1] * invWt;
        rgb1Squared[3 * offset + 2] += splatScale * pixel.L1square[2] * invWt;

        rgb2[3 * offset] += splatScale * pixel.L2[0] * invWt;
        rgb2[3 * offset + 1] += splatScale * pixel.L2[1] * invWt;
        rgb2[3 * offset + 2] += splatScale * pixel.L2[2] * invWt;
        rgb2Squared[3 * offset] += splatScale * pixel.L2square[0] * invWt;
        rgb2Squared[3 * offset + 1] += splatScale * pixel.L2square[1] * invWt;
        rgb2Squared[3 * offset + 2] += splatScale * pixel.L2square[2] * invWt;
        diff[3 * offset] += splatScale * pixel.D[0] * invWt;
        diff[3 * offset + 1] += splatScale * pixel.D[1] * invWt;
        diff[3 * offset + 2] += splatScale * pixel.D[2] * invWt;
        diffSquared[3 * offset] += splatScale * pixel.Dsquare[0] * invWt;
        diffSquared[3 * offset + 1] += splatScale * pixel.Dsquare[1] * invWt;
        diffSquared[3 * offset + 2] += splatScale * pixel.Dsquare[2] * invWt;

        recipPdfs[3 * offset] += splatScale * pixel.reciprocal_pdf * invWt;
        recipPdfs[3 * offset + 1] += splatScale * pixel.reciprocal_pdf * invWt;
        recipPdfs[3 * offset + 2] += splatScale * pixel.reciprocal_pdf * invWt;

        // Scale pixel value by _scale_
        rgb1[3 * offset] *= scale;
        rgb1[3 * offset + 1] *= scale;
        rgb1[3 * offset + 2] *= scale;
        rgb2[3 * offset] *= scale;
        rgb2[3 * offset + 1] *= scale;
        rgb2[3 * offset + 2] *= scale;
        diff[3 * offset] *= scale;
        diff[3 * offset + 1] *= scale;
        diff[3 * offset + 2] *= scale;
        
        rgb1Squared[3 * offset] *= scale * scale;
        rgb1Squared[3 * offset + 1] *= scale * scale;
        rgb1Squared[3 * offset + 2] *= scale * scale;
        rgb2Squared[3 * offset] *= scale * scale;
        rgb2Squared[3 * offset + 1] *= scale * scale;
        rgb2Squared[3 * offset + 2] *= scale * scale;
        diffSquared[3 * offset] *= scale * scale;
        diffSquared[3 * offset + 1] *= scale * scale;
        diffSquared[3 * offset + 2] *= scale * scale;

        avgRecipPdf += recipPdfs[3 * offset];
        
        ++offset;
    }
    avgRecipPdf /= offset;

    // Write RGB image
    LOG(INFO) << "Writing image " << filename << " with bounds " << croppedPixelBounds;
    pbrt::WriteBinary(filename + "_F.bin", &rgb1[0], croppedPixelBounds, fullResolution);
    pbrt::WriteBinary(filename + "_H.bin", &rgb2[0], croppedPixelBounds, fullResolution);
    pbrt::WriteBinary(filename + "_D.bin", &diff[0], croppedPixelBounds, fullResolution);
    pbrt::WriteBinary(filename + "_Fsquare.bin", &rgb1Squared[0], croppedPixelBounds, fullResolution);
    pbrt::WriteBinary(filename + "_Hsquare.bin", &rgb2Squared[0], croppedPixelBounds, fullResolution);
    pbrt::WriteBinary(filename + "_Dsquare.bin", &diffSquared[0], croppedPixelBounds, fullResolution);
    pbrt::WriteBinary(filename + "_rpdf.bin", &recipPdfs[0], croppedPixelBounds, fullResolution);

    // Divide reciprocal pdfs with ave value
    offset = 0;
    for (Point2i p : croppedPixelBounds) {
        recipPdfs[3 * offset] /= avgRecipPdf;
        recipPdfs[3 * offset + 1] /= avgRecipPdf;
        recipPdfs[3 * offset + 2] /= avgRecipPdf;
        offset++;
    }

    pbrt::WriteImage(filename + "_F.png", &rgb1[0], croppedPixelBounds, fullResolution);
    pbrt::WriteImage(filename + "_H.png", &rgb2[0], croppedPixelBounds, fullResolution);
    pbrt::WriteImage(filename + "_rpdf.png", &recipPdfs[0], croppedPixelBounds, fullResolution);
    
}

Film *CreateCvFilm(const ParamSet &params, std::unique_ptr<Filter> filter) {
    // Intentionally use FindOneString() rather than FindOneFilename() here
    // so that the rendered image is left in the working directory, rather
    // than the directory the scene file lives in.
    std::string filename = params.FindOneString("filename", "");
    if (PbrtOptions.imageFile != "") {
        if (filename != "") {
            Warning(
                "Output filename supplied on command line, \"%s\", ignored "
                "due to filename provided in scene description file, \"%s\".",
                PbrtOptions.imageFile.c_str(), filename.c_str());
        } else
            filename = PbrtOptions.imageFile;
    }
    if (filename == "") filename = "pbrt.exr";

    int xres = params.FindOneInt("xresolution", 1280);
    int yres = params.FindOneInt("yresolution", 720);
    if (PbrtOptions.quickRender) xres = std::max(1, xres / 4);
    if (PbrtOptions.quickRender) yres = std::max(1, yres / 4);
    Bounds2f crop(Point2f(0, 0), Point2f(1, 1));
    int cwi;
    const Float *cr = params.FindFloat("cropwindow", &cwi);
    if (cr && cwi == 4) {
        crop.pMin.x = Clamp(std::min(cr[0], cr[1]), 0.f, 1.f);
        crop.pMax.x = Clamp(std::max(cr[0], cr[1]), 0.f, 1.f);
        crop.pMin.y = Clamp(std::min(cr[2], cr[3]), 0.f, 1.f);
        crop.pMax.y = Clamp(std::max(cr[2], cr[3]), 0.f, 1.f);
    } else if (cr)
        Error("%d values supplied for \"cropwindow\". Expected 4.", cwi);

    Float scale = params.FindOneFloat("scale", 1.);
    Float diagonal = params.FindOneFloat("diagonal", 35.);
    Float maxSampleLuminance = params.FindOneFloat("maxsampleluminance",
                                                   Infinity);
    return new CvFilm(Point2i(xres, yres), crop, std::move(filter), diagonal,
                      filename, scale, maxSampleLuminance);    
}

}  // namespace pbrt
