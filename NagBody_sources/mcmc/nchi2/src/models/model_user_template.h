/*==============================================================================
	MODULE: models_user.h			[nchi2]
==============================================================================*/

// ==========================================
// Begin: USER Model

#define USERMODEL 100

// Begin:: Public interfaces:
global real rcUSER(real x, real params[]);
global real Chi2_USER(real (*ymodel)(real, real *), real params[]);
global void Model_USER_end(void);
// End :: Public interfaces
//


// End: USER Model
// ==========================================

