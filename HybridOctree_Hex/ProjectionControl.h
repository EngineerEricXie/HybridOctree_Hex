#ifndef PROJECTION_CONTROL_H
#define PROJECTION_CONTROL_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

struct ProjectionOptions {
    double targetScaledJacobian = 0.53;
    double surfaceTolerance = 1e-6; // In the normalized [0, 100] coordinate system.
    int maxIterations = 100000;
    int checkEvery = 1000;
    int stagnationChecks = 20;
    double progressTolerance = 1e-8;

    void Validate() const {
        if (!std::isfinite(targetScaledJacobian) || targetScaledJacobian <= 0 ||
            targetScaledJacobian >= 1 || !std::isfinite(surfaceTolerance) ||
            surfaceTolerance <= 0 || maxIterations <= 0 || checkEvery <= 0 ||
            stagnationChecks <= 0 || !std::isfinite(progressTolerance) || progressTolerance < 0)
            throw std::invalid_argument("Invalid projection stopping options");
    }
};

enum class ProjectionStatus { Running, Converged, IterationLimit, Stagnated, InvalidGeometry };

inline const char* ProjectionStatusName(ProjectionStatus status) {
    switch (status) {
    case ProjectionStatus::Running: return "running";
    case ProjectionStatus::Converged: return "converged";
    case ProjectionStatus::IterationLimit: return "iteration limit";
    case ProjectionStatus::Stagnated: return "stagnated";
    default: return "invalid geometry";
    }
}

struct ProjectionResult {
    ProjectionStatus status = ProjectionStatus::Running;
    int iterations = 0;
    double minScaledJacobian = -std::numeric_limits<double>::infinity();
    double maxSurfaceDistance = std::numeric_limits<double>::infinity();
};

class ProjectionControl {
public:
    explicit ProjectionControl(const ProjectionOptions& settings) : options(settings) {
        options.Validate();
    }

    ProjectionStatus Check(int iteration, int badElements, double minJacobian,
                           double maxDistance) {
        if (!std::isfinite(minJacobian) || !std::isfinite(maxDistance))
            return ProjectionStatus::InvalidGeometry;
        if (minJacobian > options.targetScaledJacobian && maxDistance < options.surfaceTolerance)
            return ProjectionStatus::Converged;
        if (iteration >= options.maxIterations) return ProjectionStatus::IterationLimit;
        const double epsilon = options.progressTolerance;
        if (badElements < bestBad || minJacobian > bestJacobian + epsilon ||
            maxDistance < bestDistance - epsilon) {
            staleChecks = 0;
            bestBad = std::min(bestBad, badElements);
            bestJacobian = std::max(bestJacobian, minJacobian);
            bestDistance = std::min(bestDistance, maxDistance);
        } else {
            ++staleChecks;
        }
        return staleChecks >= options.stagnationChecks ? ProjectionStatus::Stagnated
                                                       : ProjectionStatus::Running;
    }

private:
    ProjectionOptions options;
    int staleChecks = 0, bestBad = std::numeric_limits<int>::max();
    double bestJacobian = -std::numeric_limits<double>::infinity();
    double bestDistance = std::numeric_limits<double>::infinity();
};
#endif
