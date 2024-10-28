/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/parser.h>
#include <nori/scene.h>
#include <nori/camera.h>
#include <nori/block.h>
#include <nori/timer.h>
#include <nori/bitmap.h>
#include <nori/sampler.h>
#include <nori/integrator.h>
#include <nori/gui.h>
#include <tbb/parallel_for.h>
#include <tbb/blocked_range.h>
#include <tbb/task_scheduler_init.h>
#include <filesystem/resolver.h>
#include <thread>
#include <vector>

using namespace nori;

static int threadCount = -1;
static bool gui = true;

static void renderBlock(const Scene *scene, Sampler *sampler, ImageBlock &block, ImageBlock &fullBlock) {
    const Camera *camera = scene->getCamera();
    const Integrator *integrator = scene->getIntegrator();

    Point2i offset = block.getOffset();
    Vector2i size  = block.getSize();

    /* Clear the block contents */
    block.clear();

    /* For each pixel and pixel sample sample */
    for (int y=0; y<size.y(); ++y) {
        for (int x=0; x<size.x(); ++x) {
            // variance = fullBlock.getVariance(x + offset.x(), y + offset.y());
            // // cout << variance << "\n";
            // if (variance < 1.0f)
            //     continue;
            for (uint32_t i=0; i<sampler->getSampleCount(); ++i) {
                Point2f pixelSample = Point2f((float) (x + offset.x()), (float) (y + offset.y())) + sampler->next2D();
                Point2f apertureSample = sampler->next2D();

                /* Sample a ray from the camera */
                Ray3f ray;
                Color3f value = camera->sampleRay(ray, pixelSample, apertureSample);

                /* Compute the incident radiance */
                value *= integrator->Li(scene, sampler, ray);

                /* Store in the image block */
                block.put(pixelSample, value);
            }
        }
    }
}

static void render(Scene *scene, const std::string &filename) {
    const Camera *camera = scene->getCamera();
    Vector2i outputSize = camera->getOutputSize();
    scene->getIntegrator()->preprocess(scene);

    uint32_t minSamples = 512;
    float logMinSamples = log2(minSamples);

    int numIterations;
    if (scene->getSampler()->getSampleCount() <= minSamples)
        numIterations = 1;
    else
        numIterations = std::ceil(log2(scene->getSampler()->getSampleCount()) + 1 - logMinSamples);
    int totalSamples = 0;

    /* This goes with the code that stops at every iteration to save it */
    std::vector<std::unique_ptr<Sampler>> samplers((int) std::ceil(outputSize.x() / (float) NORI_BLOCK_SIZE) * (int) std::ceil(outputSize.y() / (float) NORI_BLOCK_SIZE));

    /* Create a block generator (i.e. a work scheduler) */
    /* This goes with the part of the code with the highest possible parallelism */
    // BlockGenerator blockGenerator(outputSize, NORI_BLOCK_SIZE, numIterations); 
    // std::vector<std::unique_ptr<Sampler>> samplers(blockGenerator.getBlocksLeft());

    for (size_t i=0; i < samplers.size(); i++) {
        /* Create a clone of the sampler per block in the grid */
        samplers[i] = scene->getSampler()->clone();
        // samplers[i].get()->setSampleCount(32);
    }

    /* Allocate memory for the entire output image and clear it */
    ImageBlock result(outputSize, camera->getReconstructionFilter());
    result.clear();

    /* Create a window that visualizes the partially rendered result */
    NoriScreen *screen = nullptr;
    if (gui) {
        nanogui::init();
        screen = new NoriScreen(result);
    }

    /* Do the following in parallel and asynchronously */
    std::thread render_thread([&] {
        tbb::task_scheduler_init init(threadCount);

        // std::cout << "Rendering .. \n";
        // std::cout.flush();
        Timer timer;



        /* This code renderes the image incrementally but stops at every iteration to save it */
        for (int iteration = 0; iteration < numIterations; iteration++) { 
            
            BlockGenerator blockGenerator = BlockGenerator(outputSize, NORI_BLOCK_SIZE);
            int sampleCount;
            if (iteration + 1 == 1)
                sampleCount = numIterations>1 ? minSamples : scene->getSampler()->getSampleCount();
            else if (iteration + 1 == numIterations)
                sampleCount = scene->getSampler()->getSampleCount() - totalSamples;
            else
                sampleCount = pow(2,floor(iteration+logMinSamples-1));

            totalSamples += sampleCount;

            std::cout << "Rendering " << sampleCount << " samples, " << totalSamples << "/" << scene->getSampler()->getSampleCount() << " overall\n";
            std::cout.flush();

            tbb::blocked_range<int> range(0, blockGenerator.getTotalBlocksLeft());
            Sampler* sampler;

            auto map = [&](const tbb::blocked_range<int> &range) {
                /* Allocate memory for a small image block to be rendered
                by the current thread */
                ImageBlock block(Vector2i(NORI_BLOCK_SIZE),
                    camera->getReconstructionFilter());

                for (int i=range.begin(); i<range.end(); ++i) {
                    /* Request an image block from the block generator */
                    std::pair<int,int> generatorResult = blockGenerator.next(block);

                    // float variance = result.getBlockVariance(block);
                    // if (variance < 0.00f) {
                    //     // cout << "Block discarded\n";
                    //     continue;
                    // }

                    sampler = samplers[generatorResult.first].get();
                    /* Inform the sampler about the block to be rendered if it's the first execution on this block */
                    if (iteration + 1 == 1)
                        sampler->prepare(block);
                    sampler->setSampleCount(sampleCount);

                    /* Render all contained pixels */
                    renderBlock(scene, sampler, block, result);

                    /* The image block has been processed. Now add it to
                    the "big" block that represents the entire image */
                    result.put(block);
                }
            };

            /// Default: parallel rendering
            tbb::parallel_for(range, map);

            /* Now turn the rendered image block into
            a properly normalized bitmap */
            std::unique_ptr<Bitmap> bitmap(result.toBitmap());

            /* Determine the filename of the output bitmap */
            std::cout << "\n";
            std::string outputName = filename;
            size_t lastdot = outputName.find_last_of(".");
            if (lastdot != std::string::npos)
                outputName.erase(lastdot, std::string::npos);
            outputName.append(tfm::format("_%dsamples", totalSamples));
            /* If the file exists, ask if the user wants to overwrite it */
            struct stat buffer;   
            if (stat((outputName+".exr").c_str(), &buffer) == 0 || stat((outputName+".png").c_str(), &buffer) == 0) {
                std::cout << "An EXR or PNG file named " << outputName << " already exists, do you want to overwrite it? [y/n]: ";
                std::string line;
                std::getline(std::cin, line);
                if (line.at(0) == 'y') {
                    /* Save using the OpenEXR format */
                    bitmap->saveEXR(outputName);

                    /* Save tonemapped (sRGB) output using the PNG format */
                    bitmap->savePNG(outputName);
                } else 
                    std::cout << "File was not saved\n";
            } else {
                /* Save using the OpenEXR format */
                bitmap->saveEXR(outputName);

                /* Save tonemapped (sRGB) output using the PNG format */
                bitmap->savePNG(outputName);
            }

        }




        /* This code lets the image be rendered incrementally, and with the highest possible parallelism */
        // tbb::blocked_range<int> range(0, blockGenerator.getTotalBlocksLeft());
        // Sampler* sampler;

        // auto map = [&](const tbb::blocked_range<int> &range) {
        //     /* Allocate memory for a small image block to be rendered
        //        by the current thread */
        //     ImageBlock block(Vector2i(NORI_BLOCK_SIZE),
        //         camera->getReconstructionFilter());

        //     for (int i=range.begin(); i<range.end(); ++i) {
        //         /* Request an image block from the block generator */
        //         std::pair<int,int> generatorResult = blockGenerator.next(block);

        //         // float variance = result.getBlockVariance(block);
        //         // if (variance < 0.00f) {
        //         //     // cout << "Block discarded\n";
        //         //     continue;
        //         // }

        //         sampler = samplers[generatorResult.first].get();
        //         /* Inform the sampler about the block to be rendered if it's the first execution on this block */
        //         if (generatorResult.second == 1) {
        //             sampler->prepare(block);
        //             sampler->setSampleCount(32);
        //         } else
        //             sampler->setSampleCount(pow(2,generatorResult.second+3));

        //         /* Render all contained pixels */
        //         renderBlock(scene, sampler, block, result);

        //         /* The image block has been processed. Now add it to
        //            the "big" block that represents the entire image */
        //         result.put(block);
        //     }
        // };

        // /// Default: parallel rendering
        // tbb::parallel_for(range, map);


        /// (equivalent to the following single-threaded call)
        // map(range);

        // std::cout << result.maxVariance() << "\n";
        
       std::cout << "Rendering done. (took " << timer.elapsedString() << ")" << endl;
    });

    /* Enter the application main loop */
    if (gui)
        nanogui::mainloop(50.f);

    /* Shut down the user interface */
    render_thread.join();

    if (gui) {
        delete screen;
        nanogui::shutdown();
    }

    // /* Now turn the rendered image block into
    //    a properly normalized bitmap */
    // std::unique_ptr<Bitmap> bitmap(result.toBitmap());

    // /* Determine the filename of the output bitmap */
    // std::string outputName = filename;
    // size_t lastdot = outputName.find_last_of(".");
    // if (lastdot != std::string::npos)
    //     outputName.erase(lastdot, std::string::npos);

    // /* Save using the OpenEXR format */
    // bitmap->saveEXR(outputName);

    // /* Save tonemapped (sRGB) output using the PNG format */
    // bitmap->savePNG(outputName);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "Syntax: " << argv[0] << " <scene.xml> [--no-gui] [--threads N]" <<  endl;
        return -1;
    }

    std::string sceneName = "";
    std::string exrName = "";

    for (int i = 1; i < argc; ++i) {
        std::string token(argv[i]);
        if (token == "-t" || token == "--threads") {
            if (i+1 >= argc) {
                cerr << "\"--threads\" argument expects a positive integer following it." << endl;
                return -1;
            }
            threadCount = atoi(argv[i+1]);
            i++;
            if (threadCount <= 0) {
                cerr << "\"--threads\" argument expects a positive integer following it." << endl;
                return -1;
            }

            continue;
        }
        else if (token == "--no-gui") {
            gui = false;
            continue;
        }

        filesystem::path path(argv[i]);

        try {
            if (path.extension() == "xml") {
                sceneName = argv[i];

                /* Add the parent directory of the scene file to the
                   file resolver. That way, the XML file can reference
                   resources (OBJ files, textures) using relative paths */
                getFileResolver()->prepend(path.parent_path());
            } else if (path.extension() == "exr") {
                /* Alternatively, provide a basic OpenEXR image viewer */
                exrName = argv[i];
            } else {
                cerr << "Fatal error: unknown file \"" << argv[i]
                     << "\", expected an extension of type .xml or .exr" << endl;
            }
        } catch (const std::exception &e) {
            cerr << "Fatal error: " << e.what() << endl;
            return -1;
        }
    }

    if (exrName !="" && sceneName !="") {
        cerr << "Both .xml and .exr files were provided. Please only provide one of them." << endl;
        return -1;
    }
    else if (exrName == "" && sceneName == "") {
        cerr << "Please provide the path to a .xml (or .exr) file." << endl;
        return -1;
    }
    else if (exrName != "") {
        if (!gui) {
            cerr << "Flag --no-gui was set. Please remove it to display the EXR file." << endl;
            return -1;
        }
        try {
            Bitmap bitmap(exrName);
            ImageBlock block(Vector2i((int) bitmap.cols(), (int) bitmap.rows()), nullptr);
            block.fromBitmap(bitmap);
            nanogui::init();
            NoriScreen *screen = new NoriScreen(block);
            nanogui::mainloop(50.f);
            delete screen;
            nanogui::shutdown();
        } catch (const std::exception &e) {
            cerr << e.what() << endl;
            return -1;
        }
    }
    else { // sceneName != ""
        if (threadCount < 0) {
            threadCount = tbb::task_scheduler_init::automatic;
        }
        try {
            std::unique_ptr<NoriObject> root(loadFromXML(sceneName));
            /* When the XML root object is a scene, start rendering it .. */
            if (root->getClassType() == NoriObject::EScene)
                render(static_cast<Scene *>(root.get()), sceneName);
        } catch (const std::exception &e) {
            cerr << e.what() << endl;
            return -1;
        }
    }

    return 0;
}
