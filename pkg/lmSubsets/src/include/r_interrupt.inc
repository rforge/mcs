#define MCS_CHECK_INTERRUPT()  r_interrupt_check()



const int r_interrupt_period_ = 8192;



int r_interrupt_tick_;
bool r_interrupt_flag_;



void r_interrupt_setup();
bool r_interrupt_check();
bool r_interrupt_flag();
void r_interrupt_check_aux(void *dummy);



inline void
r_interrupt_setup()
{
  r_interrupt_tick_ = r_interrupt_period_;
  r_interrupt_flag_ = false;
}



inline bool
r_interrupt_check()
{
  --r_interrupt_tick_;

  if (r_interrupt_tick_ == 0) {
    r_interrupt_tick_ = r_interrupt_period_;
    r_interrupt_flag_ = (R_ToplevelExec(r_interrupt_check_aux, NULL) == FALSE);
  }

  return r_interrupt_flag_;
}



inline bool
r_interrupt_flag()
{
  return r_interrupt_flag_;
}



inline void
r_interrupt_check_aux(void*)
{
  R_CheckUserInterrupt();
}
