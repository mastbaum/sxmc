#include <gtest/gtest.h>
#include "pdfz.h"

TEST(PdfzError, Constructor)
{
	pdfz::Error err("test");
	EXPECT_EQ(err.msg, "test");
}

TEST(SystematicObjects, ShiftSystematicConstructor)
{
	pdfz::ShiftSystematic syst(1, 0);
	EXPECT_EQ(pdfz::Systematic::SHIFT, syst.type);
	EXPECT_EQ(1, syst.obs);
	EXPECT_EQ(0, syst.par);
}

TEST(SystematicObjects, ScaleSystematicConstructor)
{
	pdfz::ScaleSystematic syst(0, 3);
	EXPECT_EQ(pdfz::Systematic::SCALE, syst.type);
	EXPECT_EQ(0, syst.obs);
	EXPECT_EQ(3, syst.par);
}

TEST(SystematicObjects, ResolutionSystematicConstructor)
{
	pdfz::ResolutionSystematic syst(0, 2, 1);
	EXPECT_EQ(pdfz::Systematic::RESOLUTION, syst.type);
	EXPECT_EQ(0, syst.obs);
	EXPECT_EQ(2, syst.true_obs);
	EXPECT_EQ(1, syst.par);
}
